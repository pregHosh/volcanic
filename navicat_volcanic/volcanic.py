#!/usr/bin/env python

import hashlib
import json
import math
import os
import sys
import traceback

import numpy as np
import pandas as pd

from navicat_volcanic import __version__, dv1, dv2
from navicat_volcanic.dv1 import curate_d, find_1_dv
from navicat_volcanic.dv2 import find_2_dv
from navicat_volcanic.exceptions import (
    ChoiceError,
    DescriptorChoiceError,
    InputError,
    MissingDataError,
)
from navicat_volcanic.helpers import (
    RUN_MODE_FLAGS,
    RUN_MODES,
    _headless,
    _strip_descriptor,
    arraydump,
    candidates_1_dv,
    candidates_2_dv,
    choose_1_dv,
    choose_2_dv,
    group_data_points,
    merge_descriptor_frames,
    outname,
    outpath,
    output_files,
    parse_args,
    set_outdir,
    setflags,
    user_choose_1_dv,
    user_choose_2_dv,
)
from navicat_volcanic.profile import role_of, validate_profile
from navicat_volcanic.plotting2d import (
    plot_2d_es_volcano,
    plot_2d_k_volcano,
    plot_2d_lsfer,
    plot_2d_t_volcano,
    plot_2d_tof_volcano,
)
from navicat_volcanic.plotting3d import (
    plot_3d_es_volcano,
    plot_3d_k_volcano,
    plot_3d_lsfer,
    plot_3d_t_volcano,
    plot_3d_tof_volcano,
)
from navicat_volcanic.tof import calc_es, calc_s_es, calc_tof

RESULT_SCHEMA = "volcanic_result/v1"
RESULT_JSON = "volcanic_result.json"
NSTDS = 3  # the outlier filter width volcanic has always used
# volcano kinds as named in volcanic_result.json: file prefix, y label, units
VOLCANO_KINDS = {
    "thermodynamic": ("t", "-ΔG(pds)", "kcal/mol"),
    "kinetic": ("k", "-ΔG(kds)", "kcal/mol"),
    "energy_span": ("es", "-δE (energy span)", "kcal/mol"),
    "tof": ("tof", "log10(TOF)", "log10(1/s)"),
}
# row-curation issues of the contract check that volcanic's own curation (outlier filter at NSTDS) replaces
_ROW_ISSUES = ("rows_dropped", "cells_imputed", "outliers", "few_rows")


def run_volcanic(arguments=None):
    """Command line entry point (python -m navicat_volcanic, navicat_volcanic).

    Returns the exit status: 0 the run finished; 1 the input was refused (a table that breaks the reaction_profile/v1
    contract, a bad descriptor file, no row left); 2 a choice is missing or invalid (no --descriptor, or no run mode,
    while stdin is not a terminal; an unknown descriptor column); 3 an internal error. Every run writes
    volcanic_result.json (in --outdir, else the working directory), also when it stops; argparse rejections of
    unknown options (status 2) happen before a run starts and write nothing.
    """
    args = parse_args(sys.argv[1:] if arguments is None else arguments)
    res = _new_result(args)
    code = 3
    try:
        set_outdir(args.outdir)
        _run(args, res)
        code = 0
    except ChoiceError as e:
        code = _fail(res, "refused", e, 2)
    except (InputError, MissingDataError) as e:
        code = _fail(res, "refused", e, 1)
    except Exception as e:
        res["traceback"] = traceback.format_exc()
        code = _fail(res, "error", e, 3)
    finally:
        res["files"] = output_files()
        _write_result(args, res)
    return code


def _new_result(args):
    return {
        "schema": RESULT_SCHEMA,
        "version": __version__,
        "status": "error",
        "error": "the run did not finish",
        "inputs": {
            "files": [],
            "T_K": args.temp,
            "nd": args.nd,
            "runmode": args.runmode,
            "runmode_name": RUN_MODES.get(args.runmode, "ask"),
            "imputer": args.imputer_strat,
            "nstds": NSTDS,
            "seed": args.seed,
            "descriptor": args.descriptor,
            "roles": args.roles,
            "npoints": args.npoints,
        },
        "validation": {"clear": None, "roles_source": None, "issues": []},
        "role_map": {},
        "rows": {"n_in": 0, "n_used": 0, "dropped": [], "imputed": []},
        "descriptor": {"chosen": [], "candidates": [], "how": ""},
        "lsr": [],
        "volcanoes": {},
        "points": {},
        "files": [],
        "warnings": [],
    }


def _fail(res, status, err, code):
    res["status"] = status
    res["error"] = f"{type(err).__name__}: {err}" if status == "error" else str(err)
    cands = getattr(err, "candidates", None)
    if isinstance(err, DescriptorChoiceError) and cands:
        if isinstance(cands[0], dict):
            res["descriptor"]["candidates"] = cands
        else:
            res["descriptor"]["candidates"] = [
                {"cols": list(c) if isinstance(c, (list, tuple)) else [c], "r2": r} for c, r in cands
            ]
    print(f"volcanic: {status}: {res['error']}", file=sys.stderr)
    if isinstance(err, DescriptorChoiceError) and res["descriptor"]["candidates"]:
        print("Candidate descriptors (best mean r2 first; see volcanic_result.json for all):", file=sys.stderr)
        for c in res["descriptor"]["candidates"][:10]:
            print(f"  {' + '.join(c['cols'])}  r2 {c.get('r2')}  MAE {c.get('mae')}", file=sys.stderr)
    return code


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return None if math.isnan(f) or math.isinf(f) else f
    return obj


def _write_result(args, res):
    path = os.path.join(args.outdir, RESULT_JSON) if args.outdir else RESULT_JSON
    try:
        if args.outdir:
            os.makedirs(args.outdir, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(_jsonable(res), fh, indent=1, ensure_ascii=False)
            fh.write("\n")
    except OSError as e:
        print(f"volcanic: could not write {path}: {e}", file=sys.stderr)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _read_frame(path):
    ext = str(path).rsplit(".", 1)[-1].lower() if "." in str(path) else ""
    if ext in ("xls", "xlsx"):
        df = pd.read_excel(path)
    elif ext == "csv":
        df = pd.read_csv(path, encoding="utf-8-sig")
    else:
        raise InputError(
            f"File termination for filename {path} was not understood. Try csv or one of ['xls', 'xlsx']."
        )
    df.columns = [str(c).strip() for c in df.columns]
    return df.dropna(how="all").reset_index(drop=True)


def _load(args, res):
    """Read and check the reaction tables (-i) with the reaction_profile/v1 contract, then append the descriptor
    files (-df). Returns (table, role map)."""
    frames, role_map = [], None
    for f in args.filenames:
        if f is None:
            raise InputError("-i needs a file name")
        if not os.path.isfile(f):
            raise InputError(f"input file not found: {f}")
        res["inputs"]["files"].append({"path": os.path.abspath(f), "sha256": _sha256(f), "kind": "reaction_profiles"})
        rep = validate_profile(f, consumer="volcanic", imputer=args.imputer_strat, seed=args.seed, roles=args.roles)
        res["validation"]["roles_source"] = rep["roles_source"]
        for issue in rep["issues"]:
            if issue["code"] in _ROW_ISSUES:
                continue
            res["validation"]["issues"].append(dict(issue, table=f))
            if issue["level"] == "warning":
                res["warnings"].append(f"{f}: {issue['message']}")
        if not rep["clear"]:
            res["validation"]["clear"] = False
            res["role_map"] = rep["role_map"]
            errors = [i["message"] for i in rep["issues"] if i["level"] == "error"]
            raise InputError(
                f"{f} was refused by the reaction_profile/v1 contract: " + "; ".join(errors)
                + " (python -m navicat_volcanic.profile validate TABLE prints the full report)"
            )
        df = _read_frame(f)
        cols = list(df.columns[1:])
        if cols != list(rep["role_map"]):
            raise InputError(f"{f}: columns read {cols} differ from the checked columns {list(rep['role_map'])}")
        if role_map is None:
            role_map = dict(rep["role_map"])
        elif role_map != rep["role_map"]:
            raise InputError(f"{f} does not have the same columns and roles as {args.filenames[0]}")
        frames.append(df)
    res["validation"]["clear"] = True
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    if args.verb > 1:
        print("Final reaction profile database (10 top rows):")
        print(df.head(10))
    ids = [str(v).strip() for v in df[df.columns[0]].tolist()]
    for f in args.dfilenames:
        if not os.path.isfile(f):
            raise InputError(f"descriptor file not found: {f}")
        res["inputs"]["files"].append({"path": os.path.abspath(f), "sha256": _sha256(f), "kind": "descriptors"})
        ddf = _read_frame(f)
        bad = [c for c in ddf.columns if not c or c.startswith("Unnamed:")]
        if bad or len(set(ddf.columns)) != len(ddf.columns):
            raise InputError(f"descriptor file {f} has empty or repeated headers: {list(ddf.columns)}")
        first = ddf.iloc[:, 0]
        if ddf.shape[1] > 1 and pd.to_numeric(first, errors="coerce").isna().any():
            # a first column of names: rows are matched to the reaction table by id, not by position
            dids = [str(v).strip() for v in first.tolist()]
            if sorted(dids) != sorted(ids) or len(set(dids)) != len(dids):
                raise InputError(f"descriptor file {f}: its first column does not hold the ids of the reaction table "
                                 "once each")
            ddf = ddf.iloc[[dids.index(i) for i in ids], 1:].reset_index(drop=True)
        elif ddf.shape[0] != df.shape[0]:
            raise InputError(f"descriptor file {f} has {ddf.shape[0]} rows and the reaction table {df.shape[0]}; "
                             "rows are read in the table's order (or add an id column)")
        for c in ddf.columns:
            vals = pd.to_numeric(ddf[c], errors="coerce")
            if (vals.isna() & ddf[c].notna()).any():
                raise InputError(f"descriptor file {f}: column {c!r} has non-numeric cells")
            ddf[c] = vals
        if args.verb > 1:
            print("Final descriptor database (top rows):")
            print(ddf.head())
        before = list(df.columns)
        df = merge_descriptor_frames(df, [ddf])
        for c in df.columns[len(before):]:
            role_map[c] = "DESCRIPTOR"
    return df, role_map


def _lsr(d, tags, regress, idxs):
    """Linear scaling relations of every regressed, non-constant column on the descriptor column(s) idxs, fitted
    over the rows used (float32 values as volcanic reads them, fitted in float64)."""
    out = []
    X = d[:, list(idxs)].astype(np.float64)
    for j in np.where(regress)[0]:
        y = d[:, j].astype(np.float64)
        if np.isclose(d[:, j].std(), 0):
            continue
        ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
        if ok.sum() < len(idxs) + 1:
            continue
        x, yy = X[ok], y[ok]
        if len(idxs) == 1:
            slope, intercept = np.polyfit(x[:, 0], yy, 1)
            pred = slope * x[:, 0] + intercept
            slope = float(slope)
        else:
            A = np.column_stack([x, np.ones(len(yy))])
            coef = np.linalg.lstsq(A, yy, rcond=None)[0]
            slope, intercept = [float(c) for c in coef[:-1]], coef[-1]
            pred = A @ coef
        sst = float(((yy - yy.mean()) ** 2).sum())
        r2 = 1.0 - float(((yy - pred) ** 2).sum()) / sst if sst > 0 else float("nan")
        out.append({"target": str(tags[j]), "slope": slope, "intercept": float(intercept), "r2": r2,
                    "n": int(ok.sum())})
    return out


def _profile_points(kind, d, tags, coeff, regress, dgr, T):
    """Per row, the volcano value of the row's own profile and its determining states [lower, upper], computed as
    the 2D plotting functions compute their points (energy span and TOF: TDI and TDTS; kinetic and thermodynamic:
    the two states of the determining step). dgr[i] is the reaction energy of row i."""
    tr, dr, cr = tags[regress], d[:, regress], coeff[regress]
    if kind == "thermodynamic":
        tr, dr = tr[~cr], dr[:, ~cr]
    ys, det = [], []
    for i in range(d.shape[0]):
        if kind == "tof":
            y = np.log10(calc_tof(dr[i, :], dgr[i], T, cr, exact=True)[0])
            _, imax, imin, _ = calc_es(dr[i, :], dgr[i], esp=True)
        elif kind == "energy_span":
            y, imax, imin, _ = calc_es(dr[i, :-1], dgr[i], esp=True)
        else:
            y, imax, imin, _ = calc_s_es(dr[i, :-1], dgr[i], esp=True)
        ys.append(float(np.float32(y)))
        det.append([str(tr[imin]), str(tr[imax])])
    return ys, det


def _record_2d(res, kind, tag, ids, px, py, det, rid=None, rb=None):
    prefix, ylabel, units = VOLCANO_KINDS[kind]
    entry = {
        "csv": outname(f"{prefix}_volcano_{tag}.csv"),
        "png": outname(f"{prefix}_volcano_{tag}.png"),
        "points_csv": outname(f"points_{prefix}_volcano_{tag}.csv"),
        "x_label": f"{tag}",
        "x_units": "kcal/mol",
        "y_label": ylabel,
        "units": units,
        "csv_columns": ["descriptor", "y", "ci95"],
    }
    if rid is not None and rb is not None:
        entry["regions"] = [{"label": str(r), "from_x": float(b)} for r, b in zip(rid, rb)]
    res["volcanoes"][kind] = entry
    res["points"][kind] = [
        {"id": i, "x": float(x), "y": float(y), "determining": dt} for i, x, y, dt in zip(ids, px, py, det)
    ]


def _record_3d(res, kind, tag1, tag2, ids, px, py, ys, det):
    prefix, ylabel, units = VOLCANO_KINDS[kind]
    entry = {
        "csv": outname(f"{prefix}_volcano_{tag1}_{tag2}.csv"),
        "png": outname(f"{prefix}_volcano_{tag1}_{tag2}.png"),
        "x_label": [f"{tag1}", f"{tag2}"],
        "x_units": "kcal/mol",
        "y_label": ylabel,
        "units": units,
        "csv_columns": ["descriptor_1", "descriptor_2", "y"],
    }
    if kind == "energy_span":
        entry["regions_csv"] = outname(f"regions_es_volcano_{tag1}_{tag2}.csv")
        entry["regions_png"] = outname(f"regions_es_volcano_{tag1}_{tag2}.png")
    res["volcanoes"][kind] = entry
    res["points"][kind] = [
        {"id": i, "x": [float(a), float(b)], "y": y, "determining": dt}
        for i, a, b, y, dt in zip(ids, px, py, ys, det)
    ]


def _descriptor_choice(args):
    """None (ask, or stop when headless), "auto", or the column name(s) given."""
    given = args.descriptor
    if not given:
        return None
    if len(given) == 1 and given[0].strip().lower() == "auto":
        return "auto"
    if any(g.strip().lower() == "auto" for g in given):
        raise DescriptorChoiceError("--descriptor auto cannot be combined with a column name")
    if len(given) != args.nd:
        raise DescriptorChoiceError(
            f"-nd {args.nd} needs {args.nd} --descriptor column(s) (or --descriptor auto), got {len(given)}: {given}"
        )
    return given[0] if args.nd == 1 else list(given)


def _run(args, res):
    verb = args.verb
    nd = args.nd
    runmode = args.runmode
    temp = args.temp
    imputer_strat = args.imputer_strat
    refill = args.refill
    dump = args.dump
    lmargin = args.lmargin
    rmargin = args.rmargin
    npoints = args.npoints
    plotmode = args.plotmode

    if nd not in (1, 2):
        raise InputError(f"Invalid number of descriptors in input!\n Accepted values are:\n [1, 2], not {nd}")
    if imputer_strat not in ("simple", "knn", "iterative", "none"):
        raise InputError(
            f"Invalid imputer strat in input!\n Accepted values are:\n ['simple', 'knn', 'iterative', 'none']"
        )
    # A run that cannot ask must be told what to build and which descriptor to use; stop before any work.
    choice = _descriptor_choice(args)
    if runmode not in RUN_MODES and _headless():
        raise ChoiceError(
            "no run mode was given and stdin is not a terminal, so volcanic cannot ask which volcano plots to "
            f"build; pass one of {RUN_MODE_FLAGS}"
        )
    if args.seed is not None:
        dv1.rng = np.random.default_rng(args.seed)
        dv2.rng = np.random.default_rng(args.seed)

    df, role_map = _load(args, res)
    res["role_map"] = role_map

    # Fill in reaction profile names/IDs from input data.
    if verb > 0:
        print(
            f"volcanic will assume that {df.columns[0]} contains names/IDs of reaction profiles."
        )
    names = df[df.columns[0]].values
    row_ids = [str(n).strip() for n in names]

    # Atttempts to group data points based on shared characters in names.
    cb, ms = group_data_points(args.ic, args.fc, names)

    # Expecting a full reaction profile with corresponding column names. Descriptors will be identified.
    tags = np.array([str(tag).strip() for tag in df.columns[1:]], dtype=object)
    d = np.float32(df.to_numpy()[:, 1:])

    # Column roles come from the reaction_profile/v1 check (a roles file, else the name rule role_of shared with
    # mikimo). The check has already refused tables without exactly one END column as the last energy column, or
    # with a reference state that is not exactly 0.0.
    roles = [role_map.get(tag) or role_of(tag) for tag in tags]

    # TS or intermediate are interpreted from column names. Coeffs is a boolean array.
    coeff = np.zeros(len(tags), dtype=bool)
    regress = np.zeros(len(tags), dtype=bool)
    for i, tag in enumerate(tags):
        if roles[i] == "TS":
            if verb > 0:
                print(f"Assuming field {tag} corresponds to a TS.")
            coeff[i] = True
            regress[i] = True
        elif roles[i] == "DESCRIPTOR":
            if verb > 0:
                print(
                    f"Assuming field {tag} corresponds to a non-energy descriptor variable."
                )
            tags[i] = _strip_descriptor(tag) or tag
            coeff[i] = False
            regress[i] = False
        elif roles[i] == "END":
            if verb > 0:
                print(f"Assuming ΔG of the reaction(s) are given in field {tag}.")
            dgr = d[:, i]
            coeff[i] = False
            regress[i] = True
        else:
            if verb > 0:
                print(f"Assuming field {tag} corresponds to a non-TS stationary point.")
            coeff[i] = False
            regress[i] = True

    # Your data might contain outliers (human error, computation error) or missing points.
    # We will attempt to curate your data automatically.
    # Dropped rows and filled cells are always reported (reaction_profile/v1).
    n_in = d.shape[0]
    before = np.copy(d)
    d_in = d  # curate_d writes the curated values into this array before it drops rows
    d, kept, ms = curate_d(
        d, regress, np.arange(n_in), ms, tags, imputer_strat, nstds=NSTDS, verb=verb, seed=args.seed
    )
    cb = cb[kept]
    rtags = tags[regress]
    dropped = []
    for i in sorted(set(range(n_in)) - set(kept.tolist())):
        missing = [str(rtags[k]) for k in np.where(np.isnan(before[i, regress]))[0]]
        outl = [k for k in np.where(~np.isnan(before[i, regress]) & np.isnan(d_in[i, regress]))[0]]
        if imputer_strat != "none":
            reason = "not imputable"
        else:
            reason = "outlier" if outl else "incomplete"
        entry = {"id": row_ids[i], "reason": reason, "missing": missing}
        if outl:
            entry["outliers"] = [{"column": str(rtags[k]), "value": float(before[i, regress][k])} for k in outl]
        dropped.append(entry)
    if dropped:
        print(
            f"{len(dropped)} of {n_in} reaction profiles were dropped (missing values, or outliers beyond {NSTDS} standard deviations, with imputer {imputer_strat}): {', '.join(e['id'] for e in dropped)}"
        )
    old, new = before[kept][:, regress], d[:, regress]
    filled = np.argwhere(np.isnan(old) | ~np.isclose(old, new, equal_nan=True))
    imputed = [
        dict(
            {"id": row_ids[kept[r]], "column": str(rtags[c]), "value": float(new[r, c])},
            **({"reason": "missing"} if np.isnan(old[r, c]) else {"reason": "outlier", "original": float(old[r, c])}),
        )
        for r, c in filled
    ]
    if len(filled):
        print(
            f"{len(filled)} cells were filled by the {imputer_strat} imputer (missing values or outliers): "
            + ", ".join(f"{names[kept[r]]}/{tags[regress][c]}" for r, c in filled)
        )
    res["rows"] = {"n_in": int(n_in), "n_used": int(len(kept)), "dropped": dropped, "imputed": imputed}
    if not len(kept):
        raise InputError(
            f"no reaction profile is left after curation: all {n_in} rows were dropped (missing values or outliers)"
        )
    ids_used = [row_ids[k] for k in kept]
    # The plotting functions read the reaction energy of curated row i as dgr[i], so dgr must hold the rows that were
    # kept (the table's column was read before curation; after a dropped row it would belong to another row).
    dgr = dgr[kept]

    if verb > 0:
        print(f"Reaction profile is given by stationary points:\n {tags[regress]}")

    # Runmode used to set up flags for volcano generation.
    t_volcano, k_volcano, es_volcano, tof_volcano = setflags(runmode)
    interactive = choice is None and not _headless()

    if nd == 1:
        # volcanic will find best non-TS descriptor variable
        if interactive:
            dvs, r2s = find_1_dv(d, tags, coeff, regress, verb)
            idx = user_choose_1_dv(dvs, r2s, tags)
            res["descriptor"]["candidates"] = candidates_1_dv(d, tags, coeff, regress, verb=-1)
            res["descriptor"]["how"] = "chosen at the prompt"
        else:
            pick = choose_1_dv(d, tags, coeff, regress, choice=choice, seed=args.seed, verb=verb)
            idx = pick["idx"]
            res["descriptor"]["candidates"] = pick["candidates"]
            res["descriptor"]["how"] = pick["how"]
            res["warnings"] += pick["warnings"]
            for w in pick["warnings"]:
                print(f"Warning: {w}")
        if verb > 0:
            print(f"Reaction profile is given by stationary points:\n {tags[regress]}")
        if idx is None:
            res["status"], res["error"] = "declined", ""
            print("No descriptor was chosen: nothing was built.")
            return
        idx = int(np.ravel(idx)[0])
        res["descriptor"]["chosen"] = [str(tags[idx])]
        print(f"Descriptor: {tags[idx]} ({res['descriptor']['how']})")
        res["lsr"] = _lsr(d, tags, regress, [idx])
        if idx is not None:
            print(f"Generating LSR plots using descriptor variable {tags[idx]}")
            if refill:
                d = plot_2d_lsfer(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
            else:
                _ = plot_2d_lsfer(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )

            volcano_headers = []
            volcano_list = []
            xint = None
            tag = tags[idx]

            if t_volcano:
                print(
                    f"Generating 2D thermodynamic volcano plot using descriptor variable {tags[idx]}"
                )
                xint, ymin, px, py, xmin, xmax, rid, rb = plot_2d_t_volcano(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("Thermodynamic volcano")
                volcano_list.append(ymin)
                det = _profile_points("thermodynamic", d, tags, coeff, regress, dgr, temp)[1]
                _record_2d(res, "thermodynamic", tag, ids_used, px, py, det, rid, rb)

            if k_volcano:
                print(
                    f"Generating 2D kinetic volcano plot using descriptor variable {tags[idx]}"
                )
                xint, ymin, px, py, xmin, xmax, rid, rb = plot_2d_k_volcano(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("Kinetic volcano")
                volcano_list.append(ymin)
                det = _profile_points("kinetic", d, tags, coeff, regress, dgr, temp)[1]
                _record_2d(res, "kinetic", tag, ids_used, px, py, det, rid, rb)

            if es_volcano:
                print(
                    f"Generating 2D energy span volcano plot using descriptor variable {tags[idx]}"
                )
                xint, ymin, px, py, xmin, xmax, rid, rb = plot_2d_es_volcano(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("ES volcano")
                volcano_list.append(ymin)
                det = _profile_points("energy_span", d, tags, coeff, regress, dgr, temp)[1]
                _record_2d(res, "energy_span", tag, ids_used, px, py, det, rid, rb)

            if tof_volcano:
                print(
                    f"Generating 2D TOF volcano plot using descriptor variable {tags[idx]}"
                )
                xint, ytof, px, py, xmin, xmax = plot_2d_tof_volcano(
                    idx,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    temp,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("TOF volcano")
                volcano_list.append(ytof)
                det = _profile_points("tof", d, tags, coeff, regress, dgr, temp)[1]
                _record_2d(res, "tof", tag, ids_used, px, py, det)

            if dump and xint is not None:
                arraydump(outpath("2d_volcanos.hdf5"), xint, volcano_list, volcano_headers)

    if nd == 2:
        # volcanic will find best non-TS combination of two descriptor variables
        if interactive:
            dvs, r2s = find_2_dv(d, tags, coeff, regress, verb)
            idx1, idx2 = user_choose_2_dv(dvs, r2s, tags)
            res["descriptor"]["candidates"] = candidates_2_dv(d, tags, coeff, regress, verb=-1)
            res["descriptor"]["how"] = "chosen at the prompt"
        else:
            pick = choose_2_dv(d, tags, coeff, regress, choice=choice, seed=args.seed, verb=verb)
            idx1, idx2 = pick["idx"]
            res["descriptor"]["candidates"] = pick["candidates"]
            res["descriptor"]["how"] = pick["how"]
            res["warnings"] += pick["warnings"]
            for w in pick["warnings"]:
                print(f"Warning: {w}")
        if idx1 is None or idx2 is None:
            res["status"], res["error"] = "declined", ""
            print("No descriptor pair was chosen: nothing was built.")
            return
        idx1, idx2 = int(np.ravel(idx1)[0]), int(np.ravel(idx2)[0])
        res["descriptor"]["chosen"] = [str(tags[idx1]), str(tags[idx2])]
        print(f"Descriptors: {tags[idx1]} and {tags[idx2]} ({res['descriptor']['how']})")
        res["lsr"] = _lsr(d, tags, regress, [idx1, idx2])
        if idx1 is not None and idx2 is not None:
            print(
                f"Generating multivariate LSR plots using combination of descriptor variables {tags[idx1]} and {tags[idx2]}"
            )
            if refill:
                d = plot_3d_lsfer(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
            else:
                _ = plot_3d_lsfer(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )

            volcano_headers = []
            volcano_list = []
            x1int = None
            x2int = None
            tag1, tag2 = tags[idx1], tags[idx2]

            if t_volcano:
                print(
                    f"Generating 3D thermodynamic volcano plot using combination of descriptor variables {tags[idx1]} and {tags[idx2]}"
                )
                x1int, x2int, grid, px, py = plot_3d_t_volcano(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("Thermodynamic volcano")
                volcano_list.append(grid)
                ys, det = _profile_points("thermodynamic", d, tags, coeff, regress, dgr, temp)
                _record_3d(res, "thermodynamic", tag1, tag2, ids_used, px, py, ys, det)

            if k_volcano:
                print(
                    f"Generating 3D kinetic volcano plot using combination of descriptor variables {tags[idx1]} and {tags[idx2]}"
                )
                x1int, x2int, grid, px, py = plot_3d_k_volcano(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("Kinetic volcano")
                volcano_list.append(grid)
                ys, det = _profile_points("kinetic", d, tags, coeff, regress, dgr, temp)
                _record_3d(res, "kinetic", tag1, tag2, ids_used, px, py, ys, det)

            if es_volcano:
                print(
                    f"Generating 3D energy span volcano plot using combination of descriptor variables {tags[idx1]} and {tags[idx2]}"
                )
                x1int, x2int, grid, px, py = plot_3d_es_volcano(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("ES volcano")
                volcano_list.append(grid)
                ys, det = _profile_points("energy_span", d, tags, coeff, regress, dgr, temp)
                _record_3d(res, "energy_span", tag1, tag2, ids_used, px, py, ys, det)

            if tof_volcano:
                print(
                    f"Generating 3D TOF volcano plot using combination of descriptor variables {tags[idx1]} and {tags[idx2]}"
                )
                x1int, x2int, grid, px, py = plot_3d_tof_volcano(
                    idx1,
                    idx2,
                    d,
                    tags,
                    coeff,
                    regress,
                    dgr,
                    temp,
                    cb,
                    ms,
                    lmargin,
                    rmargin,
                    npoints,
                    plotmode,
                    verb,
                )
                volcano_headers.append("TOF volcano")
                volcano_list.append(grid)
                ys, det = _profile_points("tof", d, tags, coeff, regress, dgr, temp)
                _record_3d(res, "tof", tag1, tag2, ids_used, px, py, ys, det)

            if dump and x1int is not None and x2int is not None:
                arraydump(
                    outpath("3d_volcanos.hdf5"),
                    np.vstack([x1int, x2int]),
                    volcano_list,
                    volcano_headers,
                )

    undefined = sorted({f"{k}/{p['id']}" for k, pts in res["points"].items() for p in pts
                        if not np.all(np.isfinite(np.asarray(p["y"], dtype=float)))})
    if undefined:
        res["warnings"].append(f"{len(undefined)} points have no finite value (written as null): "
                               + ", ".join(undefined[:20]) + (" ..." if len(undefined) > 20 else ""))
    res["status"], res["error"] = "ok", ""

    def volcanic_2d(
        runmode,
        idx,
        d,
        tags,
        coeff,
        regress,
        dgr,
        temp=298.15,
        cb=None,
        ms=None,
        lmargin=20,
        rmargin=20,
        npoints=200,
        plotmode=1,
        verb=0,
    ):
        if runmode == 5:
            raise InputError(
                "In function mode, options requiring explicit CLI input are disabled."
            )
        t_volcano, k_volcano, es_volcano, tof_volcano = setflags(runmode)
        d = plot_2d_lsfer(
            idx,
            d,
            tags,
            coeff,
            regress,
            cb,
            ms,
            lmargin,
            rmargin,
            npoints,
            plotmode,
            verb,
        )
        if t_volcano:
            return plot_2d_t_volcano(
                idx,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if k_volcano:
            return plot_2d_k_volcano(
                idx,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if es_volcano:
            return plot_2d_es_volcano(
                idx,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if tof_volcano:
            return plot_2d_tof_volcano(
                idx,
                d,
                tags,
                coeff,
                regress,
                dgr,
                temp,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )

    def volcanic_3d(
        runmode,
        idx1,
        idx2,
        d,
        tags,
        coeff,
        regress,
        dgr,
        temp=298.15,
        cb=None,
        ms=None,
        lmargin=20,
        rmargin=20,
        npoints=200,
        plotmode=1,
        verb=0,
    ):
        if runmode == 5:
            raise InputError(
                "In function mode, options requiring explicit CLI input are disabled."
            )
        t_volcano, k_volcano, es_volcano, tof_volcano = setflags(runmode)
        if t_volcano:
            return plot_3d_t_volcano(
                idx1,
                idx2,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if k_volcano:
            return plot_3d_k_volcano(
                idx1,
                idx2,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if es_volcano:
            return plot_3d_es_volcano(
                idx1,
                idx2,
                d,
                tags,
                coeff,
                regress,
                dgr,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
        if tof_volcano:
            return plot_3d_tof_volcano(
                idx1,
                idx2,
                d,
                tags,
                coeff,
                regress,
                dgr,
                temp,
                cb,
                ms,
                lmargin,
                rmargin,
                npoints,
                plotmode,
                verb,
            )
