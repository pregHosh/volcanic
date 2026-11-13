#!/usr/bin/env python

import argparse
import itertools
import math
import os
import re
import sys
from itertools import cycle

import numpy as np
import pandas as pd

from navicat_volcanic.exceptions import ChoiceError, DescriptorChoiceError, InputError

RUN_MODES = {0: "lsfer", 1: "thermodynamic", 2: "kinetic", 3: "energy_span", 4: "tof", 5: "all"}
RUN_MODE_FLAGS = "-lsfer, -thermo, -kinetic, -es, -tof or -all (or -r 0..5)"

# Where volcanic writes its files. set_outdir() points it at a directory; outpath() turns a raw name (built from
# column headers) into a safe file name there and records it, so a caller can find what a run wrote.
_OUTPUT = {"dir": None, "by_raw": {}, "by_safe": {}}
_UNSAFE = set("[]'\"\\/:*?<>|`$!&;{}")


def safe_name(name):
    """A file name from a raw name built from column headers: list-literal brackets and quotes, path, shell and
    wildcard characters, whitespace and control characters become '_', runs of '_' collapse, and the name never
    starts with '-' or '.'. The extension is kept."""
    stem, ext = os.path.splitext(str(name))
    if not re.fullmatch(r"\.[A-Za-z0-9]{1,6}", ext or "."):
        stem, ext = str(name), ""
    out = "".join("_" if (ch in _UNSAFE or ch.isspace() or ord(ch) < 32) else ch for ch in stem)
    out = re.sub(r"_+", "_", out).strip("_").lstrip("-.").strip("_")
    return (out or "output") + ext


def set_outdir(path=None):
    """Write every later volcanic output into ``path`` (created when missing); None writes into the working
    directory, as volcanic always did. Also forgets the names recorded by outpath()."""
    if path is not None:
        os.makedirs(path, exist_ok=True)
    _OUTPUT.update(dir=None if path is None else str(path), by_raw={}, by_safe={})


def outpath(name):
    """The path to write an output called ``name`` to: safe_name(name) inside the output directory. Two different
    raw names that clean to the same name get numbered suffixes instead of overwriting each other."""
    raw = str(name)
    if raw in _OUTPUT["by_raw"]:
        safe = _OUTPUT["by_raw"][raw]
    else:
        safe = safe_name(raw)
        stem, ext = os.path.splitext(safe)
        k = 1
        while safe in _OUTPUT["by_safe"] and _OUTPUT["by_safe"][safe] != raw:
            k += 1
            safe = f"{stem}_{k}{ext}"
        _OUTPUT["by_raw"][raw] = safe
        _OUTPUT["by_safe"][safe] = raw
    return os.path.join(_OUTPUT["dir"], safe) if _OUTPUT["dir"] else safe


def outname(name):
    """The file name (relative to the output directory) outpath() gave the raw name ``name``, or None."""
    return _OUTPUT["by_raw"].get(str(name))


def output_files():
    """The file names (relative to the output directory) given out by outpath() since the last set_outdir()."""
    return sorted(_OUTPUT["by_safe"])


def yesno(question):
    """Simple Yes/No Function."""
    prompt = f"{question} ? (y/n): "
    ans = input(prompt).strip().lower()
    if ans not in ["y", "n"]:
        print(f"{ans} is invalid, please try again...")
        return yesno(question)
    if ans == "y":
        return True
    return False


def bround(x, base: float = 10, type=None) -> float:
    if type == "max":
        return base * np.ceil(x / base)
    elif type == "min":
        return base * np.floor(x / base)
    else:
        tick = base * np.round(x / base)
        return tick


def group_data_points(bc, ec, names):
    try:
        groups = np.array([str(i)[bc:ec].upper() for i in names], dtype=object)
    except Exception as m:
        raise InputError(
            f"Grouping by name characters did not work. Error message was:\n {m}"
        )
    type_tags = np.unique(groups)
    cycol = cycle("bgrcmky")
    cymar = cycle("^ospXDvH")
    cdict = dict(zip(type_tags, cycol))
    mdict = dict(zip(type_tags, cymar))
    cb = np.array([cdict[i] for i in groups])
    ms = np.array([mdict[i] for i in groups])
    return cb, ms


def _headless():
    """True when nobody can answer a prompt (stdin missing or not a terminal)."""
    try:
        return sys.stdin is None or not sys.stdin.isatty()
    except (AttributeError, ValueError):
        return True


def _strip_descriptor(name):
    """A header without the word "descriptor" (any case), trimmed, as volcanic names descriptor columns."""
    s = str(name).strip()
    k = s.upper().find("DESCRIPTOR")
    return (s[:k] + s[k + 10:]).strip() if k >= 0 else s


def _tag_index(name, tags):
    names = [str(t).strip() for t in np.ravel(np.asarray(tags, dtype=object))]
    key = str(name).strip()
    if key in names:
        return names.index(key)
    bare = [_strip_descriptor(t) for t in names]
    if _strip_descriptor(key) in bare:
        return bare.index(_strip_descriptor(key))
    raise DescriptorChoiceError(f"{name!r} is not a column; columns are: {', '.join(names)}")


def _as_index(choice, tags):
    if isinstance(choice, (bool, np.bool_)):
        raise DescriptorChoiceError(f"invalid descriptor choice {choice!r}")
    if isinstance(choice, (int, np.integer)):
        if not 0 <= int(choice) < len(tags):
            raise DescriptorChoiceError(f"descriptor index {choice} is out of range 0..{len(tags) - 1}")
        return int(choice)
    if isinstance(choice, str):
        return _tag_index(choice, tags)
    raise DescriptorChoiceError(f"invalid descriptor choice {choice!r}")


def _rnum(v):
    """JSON-safe score rounded to 4 decimals, as find_1_dv / find_2_dv rank them (None for NaN)."""
    f = float(np.ravel(np.asarray(v, dtype=float))[0])
    return None if math.isnan(f) else round(f, 4)


def _sort_candidates(cands):
    return sorted(cands, key=lambda c: (c["r2"] is None, -(c["r2"] or 0.0),
                                         math.inf if c["mae"] is None else c["mae"], c["idx"]))


def _scored_1(d, tags, coeff, regress, verb):
    from navicat_volcanic.dv1 import score_1_dv

    maes, r2s, _ = score_1_dv(d, tags, coeff, regress, verb)
    names = [str(t).strip() for t in np.ravel(np.asarray(tags, dtype=object))]
    ts = np.asarray(coeff).astype(bool)
    out = [{"cols": [names[k + 1]], "idx": [k + 1], "r2": _rnum(r2s[k]), "mae": _rnum(maes[k]), "ts": bool(ts[k + 1])}
           for k in range(len(r2s))]
    return _sort_candidates(out)


def _scored_2(d, tags, coeff, regress, verb):
    from navicat_volcanic.dv2 import score_2_dv

    pairs, pcoeff, maes, r2s, _ = score_2_dv(d, tags, coeff, regress, verb)
    names = [str(t).strip() for t in np.ravel(np.asarray(tags, dtype=object))]
    out = [{"cols": [names[i + 1], names[j + 1]], "idx": [i + 1, j + 1], "r2": _rnum(r2s[k]), "mae": _rnum(maes[k]),
            "ts": bool(pcoeff[k])} for k, (i, j) in enumerate(pairs)]
    return _sort_candidates(out)


def _public(scored):
    return [{k: v for k, v in c.items() if k != "ts"} for c in scored if not c["ts"]]


def candidates_1_dv(d, tags, coeff, regress, verb=0):
    """Every non-TS column but the first scored as a single descriptor (dv1.score_1_dv), best mean r2 first:
    [{"cols": [name], "idx": [index into tags], "r2", "mae"}] (r2 and mae rounded to 4 decimals, None if undefined)."""
    return _public(_scored_1(d, tags, coeff, regress, verb))


def candidates_2_dv(d, tags, coeff, regress, verb=0):
    """Every pair of non-TS columns but the first scored as two descriptors (dv2.score_2_dv), best mean r2 first:
    [{"cols": [name1, name2], "idx": [i1, i2], "r2", "mae"}]. Like find_2_dv, it sets regress[-1] to False in place
    when the last column is constant."""
    return _public(_scored_2(d, tags, coeff, regress, verb))


def _listing(cands, n=10):
    shown = ", ".join(f"{' + '.join(c['cols'])} (r2 {c['r2']}, MAE {c['mae']})" for c in cands[:n])
    return shown + (f" and {len(cands) - n} more" if len(cands) > n else "")


def _choose(scored, tags, coeff, choice, seed, nd):
    cands = _public(scored)
    names = [str(t).strip() for t in np.ravel(np.asarray(tags, dtype=object))]
    res = {"idx": None, "chosen": [], "candidates": cands, "how": "", "tie": [], "seed": seed, "warnings": []}
    what = "descriptor" if nd == 1 else "descriptor pair"
    if choice is None:
        raise DescriptorChoiceError(
            f"a {what} must be chosen and none was given; pass 'auto' or "
            + ("a column" if nd == 1 else "two columns") + f". Candidates, best mean r2 first: {_listing(cands)}", cands)
    if isinstance(choice, str) and choice.strip().lower() == "auto":
        finite = [c for c in cands if c["r2"] is not None]
        if not finite:
            raise DescriptorChoiceError(f"no {what} candidate has a defined r2, so 'auto' cannot choose", cands)
        best = finite[0]["r2"]
        tied = sorted((c for c in finite if c["r2"] == best), key=lambda c: c["idx"])
        if len(tied) > 1:
            s = 0 if seed is None else int(seed)
            pick = tied[int(np.random.default_rng(s).integers(len(tied)))]
            res["tie"] = [" + ".join(c["cols"]) for c in tied]
            res["how"] = (f"auto: best mean r2 {best:.4f}, a tie between {', '.join(res['tie'])} (equal r2 to 4 "
                          f"decimals), broken with seed {s}")
            res["warnings"].append(f"descriptor tie: {', '.join(res['tie'])} have the same mean r2 {best:.4f}; "
                                   f"{' + '.join(pick['cols'])} was taken with seed {s}")
        else:
            pick = finite[0]
            res["how"] = f"auto: best mean r2 ({best:.4f}) among {len(cands)} candidates"
        # find_1_dv / find_2_dv also see a tie when a transition-state column scores the same; it is never a
        # candidate, but the tie is reported (the descriptor is then collinear with that TS)
        ts_tied = [" + ".join(c["cols"]) for c in scored if c["ts"] and c["r2"] == best]
        if ts_tied:
            res["tie"] = res["tie"] or [" + ".join(pick["cols"])]
            res["tie"] += ts_tied
            res["how"] += (f"; tie with transition-state column(s) {', '.join(ts_tied)} at the same r2, which are "
                           "not candidates")
            res["warnings"].append(f"descriptor tie: {' + '.join(pick['cols'])} and the transition-state column(s) "
                                   f"{', '.join(ts_tied)} have the same mean r2 {best:.4f}; transition states are "
                                   f"not candidates, so {' + '.join(pick['cols'])} was taken")
        idx = list(pick["idx"])
    else:
        given = list(choice) if isinstance(choice, (list, tuple)) else [choice]
        if len(given) != nd:
            raise DescriptorChoiceError(f"{nd} descriptor column(s) needed, got {len(given)}: {given!r}", cands)
        idx = [_as_index(g, tags) for g in given]
        if 0 in idx:
            raise DescriptorChoiceError(f"{names[0]!r} is the reference state (the first energy column) and cannot "
                                        "be a descriptor", cands)
        if len(set(idx)) != len(idx):
            raise DescriptorChoiceError(f"the two descriptors must be different columns, got {given!r}", cands)
        ts = np.asarray(coeff).astype(bool)
        for i in idx:
            if ts[i]:
                res["warnings"].append(f"descriptor {names[i]!r} is a transition-state column")
        res["how"] = "given: " + " + ".join(names[i] for i in idx)
    res["idx"] = idx[0] if nd == 1 else tuple(idx)
    res["chosen"] = [names[i] for i in idx]
    return res


def choose_1_dv(d, tags, coeff, regress, choice="auto", seed=None, verb=0):
    """Choose one descriptor column without asking anyone (the non-interactive counterpart of find_1_dv +
    user_choose_1_dv).

    d, tags, coeff, regress are what find_1_dv takes (tags[0] is the reference state). choice: "auto" takes the
    candidate with the best mean r2 (as find_1_dv ranks them, rounded to 4 decimals); a tie is broken with
    numpy.random.default_rng(seed) (seed None means 0) and reported. A column name (trimmed; "Descriptor X" and "X"
    name the same column) or an int index into tags takes that column. None raises DescriptorChoiceError listing
    the candidates. Returns {"idx": int index into tags, "chosen": [name], "candidates": candidates_1_dv(...),
    "how": str, "tie": [tied names] or [], "seed": seed, "warnings": [str]}.
    """
    return _choose(_scored_1(d, tags, coeff, regress, verb), tags, coeff, choice, seed, 1)


def choose_2_dv(d, tags, coeff, regress, choice="auto", seed=None, verb=0):
    """Choose a pair of descriptor columns without asking anyone (the counterpart of find_2_dv +
    user_choose_2_dv). As choose_1_dv; an explicit choice is a pair of names or indices into tags. "idx" is a pair
    of indices into tags (the convention of user_choose_2_dv's return value). Like find_2_dv, it sets regress[-1]
    to False in place when the last column is constant."""
    return _choose(_scored_2(d, tags, coeff, regress, verb), tags, coeff, choice, seed, 2)


def user_choose_1_dv(dvs, r2s, tags, choice=None):
    """Choose one descriptor column; returns its index into ``tags`` (or None when the user declines).

    choice: None asks on the terminal when stdin is one, and raises DescriptorChoiceError (listing the candidates)
    when it is not, so a headless run never waits on input(); "auto" takes the candidate with the best r2 (the first
    of ``dvs`` on a tie); an int is an index into ``tags``; a str is a column name of ``tags``.
    """
    if choice is not None or _headless():
        otags = np.asarray(tags, dtype=object)
        cands = [(str(np.ravel(otags[dv])[0]), float(np.ravel(r2)[0])) for dv, r2 in zip(dvs, r2s)]
        if choice is None:
            raise DescriptorChoiceError(
                "a descriptor must be chosen and stdin is not a terminal; pass a choice ('auto', a column name or an "
                "index). Candidates (column, r2): " + ", ".join(f"{c} ({r:.4f})" for c, r in cands), cands)
        if isinstance(choice, str) and choice.strip().lower() == "auto":
            if not len(dvs):
                raise DescriptorChoiceError("no descriptor candidates to choose from")
            return dvs[int(np.argmax([float(np.ravel(r)[0]) for r in r2s]))]
        return _as_index(choice, tags)
    for dv, r2 in zip(dvs, r2s):
        print(
            f"\n{tags[dv]} has been identified as a suitable descriptor variable with r2={np.round(r2,4)}."
        )
        ok = yesno("Continue using this variable")
        if ok:
            return dv
    if not ok:
        manual = yesno("Would you want to use some other descriptor variable instead")
        if manual:
            for i, tag in enumerate(tags):
                ok = yesno(f"Use {tag} as descriptor")
                if ok:
                    return i
    return None


def user_choose_2_dv(ddvs, r2s, tags, choice=None):
    """Choose a pair of descriptor columns; returns their indices into ``tags`` (or (None, None) when declined).

    choice: None asks on the terminal when stdin is one, and raises DescriptorChoiceError (listing the candidates)
    when it is not; "auto" takes the pair with the best r2 (the first of ``ddvs`` on a tie); a pair of ints is a pair
    of indices into ``tags``; a pair of str is a pair of column names of ``tags``.
    """
    if choice is not None or _headless():
        inner = [str(t) for t in np.ravel(np.asarray(tags, dtype=object))[1:]]
        cands = [((inner[dv[0]], inner[dv[1]]), float(np.ravel(r2)[0])) for dv, r2 in zip(ddvs, r2s)]
        if choice is None:
            raise DescriptorChoiceError(
                "a descriptor pair must be chosen and stdin is not a terminal; pass a choice ('auto', two column "
                "names or two indices). Candidates (columns, r2): "
                + ", ".join(f"{a} + {b} ({r:.4f})" for (a, b), r in cands), cands)
        if isinstance(choice, str) and choice.strip().lower() == "auto":
            if not len(ddvs):
                raise DescriptorChoiceError("no descriptor candidates to choose from")
            dv = ddvs[int(np.argmax([float(np.ravel(r)[0]) for r in r2s]))]
            return (dv[0] + 1, dv[1] + 1)
        if isinstance(choice, str) or len(choice) != 2:
            raise DescriptorChoiceError(f"a descriptor pair choice must be 'auto' or two columns, not {choice!r}")
        return _as_index(choice[0], tags), _as_index(choice[1], tags)
    tags = np.array(tags[1:], dtype=np.str_)
    ptags = []
    for pair in itertools.combinations(tags, r=2):
        ptags.append([pair[0], pair[1]])
    for dv, r2 in zip(ddvs, r2s):
        print(
            f"\nThe combination of {tags[dv[0]]} and {tags[dv[1]]} has been identified as a suitable combined descriptor variable with r2={np.round(r2,4)}."
        )
        ok = yesno("Continue using this combined descriptor variable")
        if ok:
            return (dv[0] + 1, dv[1] + 1)
    if not ok:
        manual = yesno(
            "Would you want to use some other descriptor variable combination instead"
        )
        if manual:
            for i, ptag in enumerate(ptags):
                ok = yesno(f"Use combination of {ptag[0]} and {ptag[1]} as descriptor")
                if ok:
                    idx1 = np.where(tags == np.str_(ptag[0]))[0][0] + 1
                    idx2 = np.where(tags == np.str_(ptag[1]))[0][0] + 1
                    return idx1, idx2
    return None, None


def build_parser():
    vbuilder = argparse.ArgumentParser(
        prog="volcanic",
        description="Build volcano plots from reaction energy profile data.",
        epilog="Remember to cite the volcanic paper: \n \nLaplaza, R., Das, S., Wodrich, M.D. et al. Constructing and interpreting volcano plots and activity maps to navigate homogeneous catalyst landscapes. Nat Protoc (2022). \nhttps://doi.org/10.1038/s41596-022-00726-2 \n \n - and enjoy!",
    )
    vbuilder.add_argument(
        "-version", "--version", action="version", version="%(prog)s 1.3.3"
    )
    runmode_arg = vbuilder.add_mutually_exclusive_group()
    vbuilder.add_argument(
        "-i",
        "--i",
        "-input",
        dest="filenames",
        nargs="?",
        action="append",
        type=str,
        required=True,
        help="Filename containing reaction energy data. See documentation for input and file formatting questions.",
    )
    vbuilder.add_argument(
        "-df",
        "--df",
        "-i2",
        "--i2",
        dest="dfilenames",
        action="append",
        type=str,
        default=[],
        help="Filename containing non-energy descriptors matching the reaction profiles. See documentation for input and file formatting questions.",
    )
    vbuilder.add_argument(
        "-nd",
        "--nd",
        dest="nd",
        type=int,
        default=1,
        help="Number of descriptor variables to use. (default: 1)",
    )
    vbuilder.add_argument(
        "-v",
        "--v",
        "--verb",
        dest="verb",
        type=int,
        default=0,
        help="Verbosity level of the code. Higher is more verbose and viceversa. Csv output files are written at every level. (default: 0)",
    )
    vbuilder.add_argument(
        "-r",
        "--r",
        dest="runmode",
        type=int,
        default=6,
        help="Defines the volcano plots to build.\n 0: LSRs\n 1: Thermodynamic\n 2: Kinetic\n 3: Energy Span\n 4: TOF\n 5: All of the above\n Other: Ask through CLI (default)",
    )
    runmode_arg.add_argument(
        "-lsfer",
        "--lsfer",
        "-lsr",
        "--lsr",
        dest="runmode",
        action="store_const",
        const=0,
        help="Set runmode to 0, building only linear scaling relationships.",
    )
    runmode_arg.add_argument(
        "-thermo",
        "--thermo",
        dest="runmode",
        action="store_const",
        const=1,
        help="Set runmode to 1, building only thermodynamic volcano.",
    )
    runmode_arg.add_argument(
        "-kinetic",
        "--kinetic",
        dest="runmode",
        action="store_const",
        const=2,
        help="Set runmode to 2, building only kinetic volcano.",
    )
    runmode_arg.add_argument(
        "-es",
        "--es",
        dest="runmode",
        action="store_const",
        const=3,
        help="Set runmode to 3, building only energy span volcano.",
    )
    runmode_arg.add_argument(
        "-tof",
        "--tof",
        dest="runmode",
        action="store_const",
        const=4,
        help="Set runmode to 4, building only TOF volcano.",
    )
    runmode_arg.add_argument(
        "-all",
        "--all",
        dest="runmode",
        action="store_const",
        const=5,
        help="Set runmode to 5, building all available volcanoes.",
    )
    vbuilder.add_argument(
        "-T",
        "-t",
        "--T",
        "--t",
        "--temp",
        "-temp",
        dest="temp",
        type=float,
        default=298.15,
        help="Temperature in K. (default: 298.15)",
    )
    vbuilder.add_argument(
        "-pm",
        "--pm",
        "-plotmode",
        "--plotmode",
        dest="plotmode",
        type=int,
        default=1,
        help="Plot mode for volcano and activity map plotting. Higher is more detailed, lower is basic. 3 includes uncertainties. (default: 1)",
    )
    vbuilder.add_argument(
        "-ic",
        "--ic",
        dest="ic",
        type=int,
        default=0,
        help="Initial character for grouping based on name for visualization. (default: 0)",
    )
    vbuilder.add_argument(
        "-fc",
        "--fc",
        dest="fc",
        type=int,
        default=2,
        help="Final character for grouping based on name for visualization. (default: 2)",
    )
    vbuilder.add_argument(
        "-rm",
        "--rm",
        dest="rmargin",
        type=int,
        default=10,
        help="Right margin to pad for visualization, in descriptor variable units. (default: 10)",
    )
    vbuilder.add_argument(
        "-lm",
        "--lm",
        dest="lmargin",
        type=int,
        default=10,
        help="Left margin to pad for visualization, in descriptor variable units. (default: 10)",
    )
    vbuilder.add_argument(
        "-np",
        "--np",
        dest="npoints",
        type=int,
        default=200,
        help="Number of grid points to use for visualization. (default: 200)",
    )
    vbuilder.add_argument(
        "-d",
        "--d",
        "-dump",
        "--dump",
        dest="dump",
        action="store_true",
        help="Flag to activate h5py dumping of data. (default: False)",
    )
    vbuilder.add_argument(
        "-is",
        "--is",
        dest="imputer_strat",
        type=str,
        default="none",
        help="Imputter to refill missing datapoints. Beta version. (default: None)",
    )
    vbuilder.add_argument(
        "-refill",
        "--refill",
        dest="refill",
        action="store_true",
        help="Refill missing values using LSRs. Beta version. (default: False)",
    )
    vbuilder.add_argument(
        "--descriptor",
        dest="descriptor",
        action="append",
        default=None,
        metavar="auto|COLUMN",
        help="Descriptor to use without asking: 'auto' (best mean r2; a tie is broken with --seed and reported) or a "
        "column header (trimmed; a 'Descriptor X' column can be named X). With -nd 2, 'auto' or the option twice. "
        "Without it volcanic asks at the terminal, and when stdin is not a terminal it stops (exit 2) listing the "
        "candidates.",
    )
    vbuilder.add_argument(
        "--outdir",
        dest="outdir",
        default=None,
        help="Directory for every file the run writes (created when missing), including volcanic_result.json. "
        "(default: the working directory)",
    )
    vbuilder.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=None,
        help="Seed for the tie-break of --descriptor auto (default 0), the iterative imputer and volcanic's other "
        "random tie-breaks (unseeded by default).",
    )
    vbuilder.add_argument(
        "--roles",
        dest="roles",
        default=None,
        help="Roles file (reaction_profile/v1) declaring each column as INT, TS, END or DESCRIPTOR. (default: "
        "<table name without its extension>.roles.json or roles.json next to the table if present, else the header name rule)",
    )
    return vbuilder


def parse_args(arguments):
    """The parsed command line (argparse Namespace); argparse exits with status 2 on an unknown option."""
    return build_parser().parse_args(arguments)


def merge_descriptor_frames(df, ddfs):
    """Append the columns of descriptor tables (one column per descriptor, rows in the profile's order, no id
    column) to the profile table, after its energy columns, as "Descriptor <header>" columns (a header that already
    holds the word "descriptor" is kept)."""
    if not ddfs:
        return df
    ddf = pd.concat(ddfs, axis=1) if len(ddfs) > 1 else ddfs[0]
    if df.shape[0] != ddf.shape[0]:
        raise InputError(
            f"Different number of entries in reaction profile input file ({df.shape[0]}) and descriptor file "
            f"({ddf.shape[0]}). Exiting."
        )
    df = df.copy()
    for column in ddf.columns:
        name = str(column).strip()
        if "descriptor" not in name.lower():
            name = f"Descriptor {name}"
        if name in [str(c).strip() for c in df.columns]:
            raise InputError(f"Descriptor column {name!r} is already in the reaction profile table. Exiting.")
        df[name] = ddf[column].values
    return df


def processargs(arguments):
    args = parse_args(arguments)

    dfs, ddfs = check_input(
        args.filenames,
        args.dfilenames,
        args.temp,
        args.nd,
        args.imputer_strat,
        args.verb,
    )
    if len(dfs) > 1:
        df = pd.concat(dfs)
    elif len(dfs) == 0:
        raise InputError("No input profiles detected in file. Exiting.")
    else:
        df = dfs[0]
    assert isinstance(df, pd.DataFrame)
    if args.verb > 1:
        print("Final reaction profile database (10 top rows):")
        print(df.head(10))

    if ddfs:
        if args.verb > 1:
            for ddf in ddfs:
                print("Final descriptor database (top rows):")
                print(ddf.head())
        df = merge_descriptor_frames(df, ddfs)
    return (
        df,
        args.nd,
        args.verb,
        args.runmode,
        args.temp,
        args.imputer_strat,
        args.refill,
        args.dump,
        args.ic,
        args.fc,
        args.lmargin,
        args.rmargin,
        args.npoints,
        args.plotmode,
    )


def check_input(filenames, dfilenames, temp, nd, imputer_strat, verb):
    accepted_excel_terms = ["xls", "xlsx"]
    accepted_imputer_strats = ["simple", "knn", "iterative", "none"]
    accepted_nds = [1, 2]
    dfs = []
    ddfs = []
    for filename in filenames:
        if filename.split(".")[-1] in accepted_excel_terms:
            dfs.append(pd.read_excel(filename))
        elif filename.split(".")[-1] == "csv":
            dfs.append(pd.read_csv(filename))
        else:
            raise InputError(
                f"File termination for filename {filename} was not understood. Try csv or one of {accepted_excel_terms}."
            )
    for dfilename in dfilenames:
        if dfilename.split(".")[-1] in accepted_excel_terms:
            ddfs.append(pd.read_excel(dfilename))
        elif dfilename.split(".")[-1] == "csv":
            ddfs.append(pd.read_csv(dfilename))
        else:
            raise InputError(
                f"File termination for filename {dfilename} was not understood. Try csv or one of {accepted_excel_terms}."
            )
    if not isinstance(temp, float):
        raise InputError("Invalid temperature input! Should be a float. Exiting.")
    if imputer_strat not in accepted_imputer_strats:
        raise InputError(
            f"Invalid imputer strat in input!\n Accepted values are:\n {accepted_imputer_strats}"
        )
    if not isinstance(verb, int):
        raise InputError("Invalid verbosity input! Should be a positive integer or 0.")
    if nd not in accepted_nds:
        raise InputError(
            f"Invalid number of descriptors in input!\n Accepted values ar:\n {nd}"
        )
    return dfs, ddfs


def arraydump(path: str, descriptor: np.array, volcano_list, volcano_headers):
    """Dump array as an hdf5 file."""
    import h5py

    h5 = h5py.File(path, "w")
    assert len(volcano_list) == len(volcano_headers)
    # hdf5 file is like a dictionary, every dataset has a key and a data value (which can be an array)
    h5.create_dataset("Descriptor", data=descriptor)
    for i, j in zip(volcano_list, volcano_headers):
        h5.create_dataset(j, data=i)
    h5.close()


def arrayread(path: str):
    """Read hdf5 dataset."""
    import h5py

    h5 = h5py.File(path, "r")
    volcano_headers = []
    volcano_list = []
    for key in h5.keys():
        if key == "Descriptor":
            descriptor = h5[key]
        else:
            volcano_headers.append(key)
            volcano_list.append(h5[key][()])
    return descriptor, volcano_list, volcano_headers


def setflags(runmode):
    if runmode == 0:
        t_volcano = False
        k_volcano = False
        es_volcano = False
        tof_volcano = False
    elif runmode == 1:
        t_volcano = True
        k_volcano = False
        es_volcano = False
        tof_volcano = False
    elif runmode == 2:
        t_volcano = False
        k_volcano = True
        es_volcano = False
        tof_volcano = False
    elif runmode == 3:
        t_volcano = False
        k_volcano = False
        es_volcano = True
        tof_volcano = False
    elif runmode == 4:
        t_volcano = False
        k_volcano = False
        es_volcano = False
        tof_volcano = True
    elif runmode == 5:
        t_volcano = True
        k_volcano = True
        es_volcano = True
        tof_volcano = True
    else:
        if _headless():
            raise ChoiceError(
                "no run mode was given and stdin is not a terminal, so volcanic cannot ask which volcano plots to "
                f"build; pass one of {RUN_MODE_FLAGS}"
            )
        t_volcano = yesno("Generate thermodynamic volcano plot")
        k_volcano = yesno("Generate kinetic volcano plot")
        es_volcano = yesno("Generate energy span volcano plot")
        tof_volcano = yesno("Generate TOF volcano plot")
    return t_volcano, k_volcano, es_volcano, tof_volcano


def test_filedump():
    dv = np.linspace(0, 50, 1000)
    tv = np.linspace(-34, 13, 1000)
    kv = np.linspace(15, 25, 1000)
    volcano_list = [tv, kv]
    volcano_headers = ["Thermodynamic volcano", "Kinetic volcano"]
    arraydump("hdf5_test.hdf5", dv, volcano_list, volcano_headers)
    dv, volcano_list, volcano_headers = arrayread("hdf5_test.hdf5")
    assert np.allclose(
        tv, volcano_list[volcano_headers.index("Thermodynamic volcano")], 4
    )
    assert np.allclose(kv, volcano_list[volcano_headers.index("Kinetic volcano")], 4)
    os.remove("hdf5_test.hdf5")


if __name__ == "__main__":
    test_filedump()
