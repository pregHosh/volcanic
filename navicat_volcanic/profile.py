#!/usr/bin/env python
"""Reaction-profile contract ``reaction_profile/v1``: one reading of a reaction table for volcanic and mikimo.

Both packages read the same tables (one row per catalyst, one column per stationary point). This module is the
single place that decides how such a table is read, which role each column has, and which rows are used. It never
fixes a table silently: every refusal and every change it applies is listed in ``issues``, and every dropped row and
imputed cell is listed in ``rows``.

Public API
----------
validate_profile(path, consumer="volcanic", imputer="none", seed=None, nstds=None, network=None, roles=None) -> dict
    path      reaction table: UTF-8 csv (a BOM is stripped) or the first sheet of an xlsx file.
    consumer  "volcanic": exactly one END column, and it is the last energy column (a multi-pathway table is refused);
              "mikimo": one END per pathway, and the last energy column is an END.
    imputer   "none" (default): every row with a missing value is dropped and listed. "knn", "simple" or "iterative"
              fill missing energies with volcanic's imputer (``dv1.call_imputer``, row by row as ``curate_d`` does);
              every filled cell is listed. Rows with a missing descriptor value are always dropped.
    seed      random state for the "iterative" imputer (knn and simple are deterministic).
    nstds     None (default): no outlier filter. A number: volcanic's outlier rule, NaN-aware, per energy column: a
              value is an outlier when |x| > |mean| + nstds*sd or |x| < |mean| - nstds*sd, with mean and sd (ddof 0)
              over the non-missing values; columns with sd ~ 0 have no outliers. An outlier cell becomes missing, so
              with imputer "none" its row is dropped (reason "outlier"), otherwise the cell is imputed (reason
              "outlier", original value kept in the entry).
    network   path of a mikimo rxn_network csv/xlsx: its species are matched to the table by exact name (see below).
    roles     None: use ``<table name without its extension>.roles.json`` or ``roles.json`` next to the table when
              one exists, else the legacy name rule; a path or a dict: that roles file; False: always the legacy
              name rule.
    returns   {"schema": "reaction_profile/v1", "clear": bool, "consumer": str, "table": str, "id_column": str,
               "issues": [{"level": "error"|"warning"|"fixed", "code": str, "message": str, ...names...}],
               "role_map": {trimmed header: "INT"|"TS"|"END"|"DESCRIPTOR"} (table order, column 0 left out),
               "roles_source": "legacy name rule" | <roles file> | "roles argument",
               "pathways": [[energy columns of pathway 1, ending at its END], ...],
               "rows": {"n_in": int, "n_used": int,
                        "dropped": [{"id", "reason": "incomplete"|"outlier"|"not imputable", "missing": [cols],
                                     "outliers": [{"column", "value"}]}],
                        "imputed": [{"id", "column", "value", "reason": "missing"|"outlier", "original"?}]},
               "imputer": {"strategy", "seed"}, "nstds": float|None,
               "network": {...} (only with ``network``, see below)}
    ``clear`` is False when any issue has level "error" or no row is left. Issues name what they are about in the
    fields "column", "columns", "id", "ids" or "species". It never raises on a bad table (an unreadable file is an
    error issue); it raises ArgumentError (a ValueError) only for a bad argument.

load_profile(...same arguments...) -> (report, frame)
    The same report, plus a pandas DataFrame of the rows used (column 0 the ids, then the role-map columns with
    trimmed headers, imputed cells filled), or None when the table could not be read. Check report["clear"] first.

validate_kinetic_data(path, network=None, imputer="none", seed=None, nstds=None, units="log10_k") -> dict
    A mikimo kinetic_data table. Columns: ``kf_<i>`` (forward) and ``kr_<i>`` (reverse), i = 1..n, holding
    log10 of the rate constants; the legacy names ``k_<i>`` and ``k_-<i>`` are read as kf_i and kr_i (reported).
    A header containing "descriptor" (any case) is a DESCRIPTOR before any rate-constant rule. Any other header is
    refused. Forward and reverse must both cover 1..n; with ``network`` the step count must equal n. The report has
    the same shape (role_map values "KF"|"KR"|"DESCRIPTOR"; no pathways) plus
    "kinetic": {"units", "n_steps", "forward": [cols for i = 1..n], "reverse": [...], "descriptors": [...],
    "legacy_names": bool}. Read rate constants by these names, never by position.

role_of(header) -> "INT" | "TS" | "END" | "DESCRIPTOR"
    The one legacy name rule shared by volcanic and mikimo, used when no roles file declares a column:
      DESCRIPTOR  the header contains "descriptor" (any case);
      TS          "TS" in any case, not preceded by a letter and not followed by a lowercase letter
                  (ts2, Ts1, TS(1,A2), dGRRS(TSOA), B-TS1, RO-TS1 are TS; OTs, adducts, Tsuji are not);
      END         the header contains "prod" (any case), or starts with "P-" or "P_", or is just "P"
                  (P2, Pd0L2 are not ENDs);
      INT         everything else.

pathways(role_map) -> [[column, ...], ...]
    The energy columns in table order, split after each END (one list per pathway).

save_roles(report_or_role_map, path)
    Write a resolved role map as a roles file: {"schema": "reaction_profile/v1", "roles": {column: role}}.

Table format
    Headers are trimmed (each trim is reported). Column 0 holds the row ids, whatever its header (blank allowed);
    ids must be present and unique. Columns with an empty header, duplicate headers (after trimming) and
    non-numeric cells are refused. Fully empty rows are skipped (reported). Energy columns are in cycle order, in
    kcal/mol. The first energy column is the reference state: it must be an INT and exactly 0.0 in every row; a table
    that breaks this is refused, never zeroed.

roles.json
    {"column header": "INT"|"TS"|"END"|"DESCRIPTOR", ...}, optionally wrapped as {"roles": {...}}. Columns it does not
    list fall back to the legacy rule (each one reported as a warning); a listed column that is not in the table, or
    an unknown role, is refused. The wrapped form can also declare the role of network species by name:
    {"roles": {...}, "species": {"L": "reactant"|"product"|"target"|"intermediate"|"spectator", ...}}; it is read
    only when a network is validated, takes precedence over the name rule below, and a species it names that is not
    in the network, or an unknown species role, is refused. In the wrapped form any section other than "roles", "species"
    and "schema" (which must name this contract, as the file written by --write-roles does), and a "species" section
    that is not an object, is refused too (a misspelt key is never ignored).

Network (mikimo)
    The rxn_network table: column 0 holds step labels, the other headers are species, cells are integer
    stoichiometric coefficients (empty = 0), and the last row is labelled "c0", "initial_conc" or "initial conc"
    (any case) with non-negative initial concentrations. Empty or duplicate species headers are refused. Each species
    gets a role, by exact name:
      a role declared in roles.json ("species") -> that role;
      an INT column of the table          -> "intermediate" ("catalyst" for species 0);
      R, R-x, R_x, R<digit>x              -> "reactant";
      "prod" in the name, P-x, P_x, P<digit>x, P* -> "product" ("target" when the name holds "*");
      a species with all-zero coefficients -> "spectator" (warning);
      anything else (a TS column, or a name that is neither an INT column nor a reactant/product, and not declared) is
      refused.
    Species 0 must be the table's first energy column (the catalyst resting state), and every INT column of the
    table must be a species of the network. Report section: {"path", "n_steps", "species": [{"name", "role",
    "c0"}], "targets": [names]}.

Related functions outside this module (navicat_volcanic.helpers), for choosing descriptors without prompts:
    choose_1_dv(d, tags, coeff, regress, choice="auto", seed=None, verb=0) and choose_2_dv(...) return
    {"idx", "chosen", "candidates", "how", "tie", "seed", "warnings"}; user_choose_1_dv / user_choose_2_dv take
    choice=... and never call input() when stdin is not a terminal. See their docstrings.

CLI
    python -m navicat_volcanic.profile validate TABLE [--consumer volcanic|mikimo] [--network F] [--imputer I]
        [--nstds X] [--seed N] [--roles F] [--kinetic] [--write-roles F] [--json]
    prints the report (one JSON object with --json, else a short summary); exit 0 when clear, 1 when not, 2 for a
    usage error, 3 for an internal error. Never prints a traceback.
"""

import argparse
import codecs
import csv
import io
import json
import math
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

SCHEMA = "reaction_profile/v1"
ROLES = ("INT", "TS", "END", "DESCRIPTOR")
ENERGY_ROLES = ("INT", "TS", "END")
SPECIES_ROLES = ("reactant", "product", "target", "intermediate", "spectator")
CONSUMERS = ("volcanic", "mikimo")
IMPUTERS = ("none", "knn", "simple", "iterative")
UNITS = ("log10_k",)
C0_LABELS = ("c0", "initial_conc", "initial conc")
LEGACY = "legacy name rule"


class ArgumentError(ValueError):
    """A bad argument value (consumer, imputer, nstds): the only error the validators raise. The command line
    exits with status 2 for it."""

_TS_RE = re.compile(r"(?<![A-Za-z])[Tt][Ss](?![a-z])")
_END_RE = re.compile(r"^[Pp](?:[-_]|$)")
_REACTANT_RE = re.compile(r"^[Rr](?:[-_\d]|$)")
_PRODUCT_RE = re.compile(r"^[Pp](?:[-_\d*]|$)")
_KF_RE = re.compile(r"^kf_?(\d+)$", re.IGNORECASE)
_KR_RE = re.compile(r"^kr_?(\d+)$", re.IGNORECASE)
_K_LEGACY_F_RE = re.compile(r"^k_?(\d+)$", re.IGNORECASE)
_K_LEGACY_R_RE = re.compile(r"^k_?-(\d+)$", re.IGNORECASE)

__all__ = ["validate_profile", "load_profile", "validate_kinetic_data", "role_of", "pathways", "save_roles", "SCHEMA"]


# ---------------------------------------------------------------------------------------------------- roles


def role_of(header):
    """The legacy name rule shared by volcanic and mikimo (see the module docstring)."""
    h = str(header).strip()
    if "descriptor" in h.lower():
        return "DESCRIPTOR"
    if _TS_RE.search(h):
        return "TS"
    if "prod" in h.lower() or _END_RE.match(h):
        return "END"
    return "INT"


def pathways(role_map):
    """Energy columns of a role map in table order, split after each END."""
    out, cur = [], []
    for col, role in role_map.items():
        if role not in ENERGY_ROLES:
            continue
        cur.append(col)
        if role == "END":
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


def save_roles(report_or_role_map, path):
    """Write a resolved role map (or a report's role map) as a roles file."""
    rm = report_or_role_map.get("role_map", report_or_role_map) if isinstance(report_or_role_map, dict) else None
    if not isinstance(rm, dict):
        raise ValueError("save_roles needs a report or a {column: role} dict")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"schema": SCHEMA, "roles": rm}, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    return str(path)


def _kinetic_role(header):
    """("DESCRIPTOR", None) | ("KF", i) | ("KR", i) | (None, None); third item True for a legacy k_i / k_-i name."""
    h = str(header).strip()
    if "descriptor" in h.lower():
        return "DESCRIPTOR", None, False
    for rx, role, legacy in ((_KF_RE, "KF", False), (_KR_RE, "KR", False), (_K_LEGACY_R_RE, "KR", True),
                             (_K_LEGACY_F_RE, "KF", True)):
        m = rx.match(h)
        if m:
            return role, int(m.group(1)), legacy
    return None, None, False


# ---------------------------------------------------------------------------------------------------- helpers


def _issue(level, code, message, **names):
    it = {"level": level, "code": code, "message": message}
    it.update({k: v for k, v in names.items() if v is not None})
    return it


def _num(v):
    """JSON-safe float (None for NaN)."""
    if v is None:
        return None
    f = float(v)
    return None if math.isnan(f) else f


def _blank(v):
    if v is None:
        return True
    if isinstance(v, float) and math.isnan(v):
        return True
    if isinstance(v, str) and not v.strip():
        return True
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def _short(items, n=12):
    items = list(items)
    return ", ".join(str(i) for i in items[:n]) + (f" (+{len(items) - n} more)" if len(items) > n else "")


def _read_table(path):
    """Read a csv (BOM stripped) or the first xlsx sheet.

    Returns (raw_headers, body, fixes, error): raw_headers are the header cells as written (str, '' for empty),
    body is the pandas DataFrame of the data rows with positional integer column labels, fixes is a list of issues,
    error is an error issue or None.
    """
    p = str(path)
    if not os.path.isfile(p):
        return None, None, [], _issue("error", "unreadable", f"file not found: {p}")
    ext = os.path.splitext(p)[1].lower()
    fixes = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            if ext in (".xlsx", ".xlsm", ".xls"):
                head = pd.read_excel(p, sheet_name=0, header=None, nrows=1)
                raw = ["" if _blank(v) else str(v) for v in head.iloc[0].tolist()] if len(head) else []
                df = pd.read_excel(p, sheet_name=0)
            elif ext in (".csv", ".txt", ""):
                data = open(p, "rb").read()
                text = data.decode("utf-8-sig")
                if data.startswith(codecs.BOM_UTF8):
                    fixes.append(_issue("fixed", "bom_stripped", "UTF-8 byte-order mark stripped"))
                raw = next(csv.reader(io.StringIO(text)), [])
                df = pd.read_csv(p, encoding="utf-8-sig")
            else:
                return None, None, [], _issue("error", "unreadable", f"unsupported file type {ext!r}: use csv or xlsx")
    except UnicodeDecodeError:
        return None, None, [], _issue("error", "unreadable", f"{p} is not UTF-8 text")
    except Exception as e:  # pandas/openpyxl/csv parse errors
        return None, None, [], _issue("error", "unreadable", f"cannot read {p}: {type(e).__name__}: {e}")
    if not isinstance(df.index, pd.RangeIndex):
        return None, None, [], _issue("error", "ragged", "data rows have more cells than the header row")
    width = df.shape[1]
    if len(raw) < width:
        raw = raw + [""] * (width - len(raw))
    elif len(raw) > width:
        if any(str(h).strip() for h in raw[width:]):
            return None, None, [], _issue("error", "ragged", "the header row has more cells than the data")
        raw = raw[:width]
    if width == 0:
        return None, None, [], _issue("error", "empty", "the table has no columns")
    df.columns = list(range(width))
    # whitespace-only text cells are empty cells
    for j in range(width):
        if df[j].dtype == object:
            df[j] = df[j].map(lambda v: np.nan if isinstance(v, str) and not v.strip() else v)
    empty = df.isna().all(axis=1).to_numpy()
    if empty.any():
        fixes.append(_issue("fixed", "empty_rows_skipped", f"{int(empty.sum())} fully empty rows skipped"))
        df = df.loc[~empty].reset_index(drop=True)
    return [str(h) for h in raw], df, fixes, None


def _headers(raw, issues, what="table"):
    """Trim headers, report trims, refuse empty and duplicate ones. Returns (trimmed, ok_columns)."""
    hdr = [h.strip() for h in raw]
    for j in range(1, len(raw)):
        if hdr[j] and raw[j] != hdr[j]:
            issues.append(_issue("fixed", "header_trimmed", f"{what} header {raw[j]!r} trimmed to {hdr[j]!r}",
                                 column=hdr[j]))
    seen = {}
    for j in range(1, len(hdr)):
        if not hdr[j]:
            issues.append(_issue("error", "empty_header",
                                 f"{what} column {j} has an empty header (an unnamed column, e.g. from trailing commas)",
                                 column=f"Unnamed: {j}"))
        else:
            seen.setdefault(hdr[j], []).append(j)
    for h, js in seen.items():
        if len(js) > 1:
            issues.append(_issue("error", "duplicate_header",
                                 f"{what} header {h!r} appears {len(js)} times (columns {js}); every header must be unique",
                                 column=h))
    ok = [j for j in range(1, len(hdr)) if hdr[j]]
    return hdr, ok


def _load_roles(roles, table_path):
    """Returns (declared {column: role} or None, source, issues, declared {species: role})."""
    if roles is False:
        return None, LEGACY, [], {}
    if roles is None:
        base = os.path.splitext(str(table_path))[0]
        for cand in (base + ".roles.json", os.path.join(os.path.dirname(str(table_path)), "roles.json")):
            if os.path.isfile(cand):
                roles = cand
                break
        else:
            return None, LEGACY, [], {}
    if isinstance(roles, dict):
        obj, source = roles, "roles argument"
    else:
        source = str(roles)
        try:
            with open(source, encoding="utf-8-sig") as fh:
                obj = json.load(fh)
        except Exception as e:
            return None, source, [_issue("error", "roles_unreadable", f"cannot read roles file {source}: {e}")], {}
    raw_species, wrap_issues = None, []
    if isinstance(obj, dict) and (isinstance(obj.get("roles"), dict) or isinstance(obj.get("species"), dict)):
        raw_species = obj.get("species")
        for key in obj:
            if key not in ("roles", "species", "schema"):
                wrap_issues.append(_issue("error", "roles_unreadable", f"roles file {source} has the unknown section "
                                          f"{key!r} (use \"roles\" and \"species\")"))
        if "schema" in obj and obj["schema"] != SCHEMA:
            wrap_issues.append(_issue("error", "roles_unreadable", f"roles file {source} declares the schema "
                                      f"{obj['schema']!r}; this reader understands {SCHEMA!r}"))
        if raw_species is not None and not isinstance(raw_species, dict):
            wrap_issues.append(_issue("error", "roles_unreadable", f"roles file {source}: \"species\" must be an "
                                      "object {species: role}"))
        obj = obj["roles"] if isinstance(obj.get("roles"), dict) else {}
    if not isinstance(obj, dict):
        return None, source, [_issue("error", "roles_unreadable",
                                     f"roles file {source} is not a {{column: role}} object")], {}
    out, issues, species = {}, list(wrap_issues), {}
    for k, v in (raw_species if isinstance(raw_species, dict) else {}).items():
        name, role = str(k).strip(), str(v).strip().lower()
        if role not in SPECIES_ROLES:
            issues.append(_issue("error", "bad_species_role", f"roles file gives network species {name!r} the unknown "
                                 f"role {v!r} (use one of {', '.join(SPECIES_ROLES)})", species=name))
            continue
        species[name] = role
    for k, v in obj.items():
        name, role = str(k).strip(), str(v).strip().upper()
        if role not in ROLES + ("ID",):
            issues.append(_issue("error", "bad_role", f"roles file gives column {name!r} the unknown role {v!r} "
                                 f"(use one of {', '.join(ROLES)})", column=name))
            continue
        out[name] = role
    if species and not out:
        out = None  # a file that only declares species leaves every column to the name rule, without a warning each
    return out, source, issues, species


def _numeric(body, cols, hdr, ids, issues):
    """Float matrix of the given columns; non-numeric cells are refused (and read as missing)."""
    out = np.full((body.shape[0], len(cols)), np.nan)
    bad = []
    for k, j in enumerate(cols):
        s = body[j]
        conv = pd.to_numeric(s, errors="coerce")
        for i in np.where(s.notna().to_numpy() & conv.isna().to_numpy())[0]:
            bad.append((ids[i], hdr[j], str(s.iloc[i])))
        out[:, k] = conv.to_numpy(dtype=float)
    for cid, col, val in bad[:20]:
        issues.append(_issue("error", "non_numeric", f"row {cid!r}, column {col!r}: {val!r} is not a number",
                             id=cid, column=col))
    if len(bad) > 20:
        issues.append(_issue("error", "non_numeric", f"{len(bad) - 20} more non-numeric cells"))
    return out


def _ids(body, issues):
    raw = body[0].tolist() if body.shape[1] else []
    ids = [str(v) for v in raw]
    missing = [i for i, v in enumerate(raw) if _blank(v)]
    if missing:
        issues.append(_issue("error", "missing_id", f"{len(missing)} rows have no id in column 0 "
                             f"(data rows {_short(i + 1 for i in missing)})"))
    seen, dup = set(), []
    for i, v in enumerate(ids):
        if i in missing:
            continue
        if v in seen and v not in dup:
            dup.append(v)
        seen.add(v)
    if dup:
        issues.append(_issue("error", "duplicate_id", f"ids in column 0 must be unique; repeated: {_short(dup)}",
                             ids=dup))
    return ids


def _outlier_cells(X, nstds):
    """{(row, col)} flagged by volcanic's rule, NaN-aware (see the module docstring)."""
    flagged = set()
    for k in range(X.shape[1]):
        col = X[:, k]
        ok = ~np.isnan(col)
        if ok.sum() < 2:
            continue
        mean = float(np.mean(col[ok]))
        sd = float(np.std(col[ok]))
        if sd <= 1e-9 * max(1.0, abs(mean)):
            continue
        moe = nstds * sd
        a = np.abs(col)
        with np.errstate(invalid="ignore"):
            hit = ok & ((a > abs(mean) + moe) | (a < abs(mean) - moe))
        for i in np.where(hit)[0]:
            flagged.add((int(i), k))
    return flagged


def _impute(X, strategy, seed):
    """volcanic's imputation, row by row as dv1.curate_d does it (rows already filled serve as donors)."""
    from navicat_volcanic.dv1 import call_imputer

    X = X.copy()
    n = X.shape[0]
    for j in range(n):
        if not np.isnan(X[j]).any():
            continue
        ref = X[np.arange(n) != j]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                new = np.asarray(call_imputer(X[j].copy(), ref, strategy, seed=seed), dtype=float).ravel()
        except Exception:
            continue
        if new.shape == X[j].shape:
            X[j] = new
    return X


def _curate(ids, V, value_k, desc_k, names, imputer, seed, nstds, issues):
    """Outlier filter, missing values and imputation. Returns (rows dict, kept mask, filled V)."""
    n = V.shape[0]
    X = V[:, value_k].copy() if value_k else np.zeros((n, 0))
    vnames = [names[k] for k in value_k]
    outl = _outlier_cells(X, nstds) if nstds is not None and X.size else set()
    orig = X.copy()
    for i, k in outl:
        X[i, k] = np.nan
    missing0 = np.isnan(orig)
    if imputer != "none" and X.size and np.isnan(X).any():
        X = _impute(X, imputer, seed)
    dropped, imputed = [], []
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        out_cols = sorted(k for (r, k) in outl if r == i)
        miss = [vnames[k] for k in range(len(vnames)) if missing0[i, k]]
        miss += [names[k] for k in desc_k if np.isnan(V[i, k])]
        left = [vnames[k] for k in range(len(vnames)) if np.isnan(X[i, k])]
        desc_gap = [names[k] for k in desc_k if np.isnan(V[i, k])]
        if left or desc_gap:
            keep[i] = False
            if imputer == "none":
                reason = "outlier" if out_cols else "incomplete"
            else:
                reason = "not imputable" if left else "incomplete"
            entry = {"id": ids[i], "reason": reason, "missing": miss}
            if out_cols:
                entry["outliers"] = [{"column": vnames[k], "value": _num(orig[i, k])} for k in out_cols]
            dropped.append(entry)
            continue
        for k in range(len(vnames)):
            if (i, k) in outl:
                imputed.append({"id": ids[i], "column": vnames[k], "value": _num(X[i, k]), "reason": "outlier",
                                "original": _num(orig[i, k])})
            elif missing0[i, k]:
                imputed.append({"id": ids[i], "column": vnames[k], "value": _num(X[i, k]), "reason": "missing"})
    Vf = V.copy()
    if value_k:
        Vf[:, value_k] = X
    if outl:
        issues.append(_issue("fixed", "outliers", f"{len(outl)} outlier cells found (nstds {nstds}): "
                             f"{_short(f'{ids[i]}/{vnames[k]}' for i, k in sorted(outl))}"))
    if dropped:
        by = {}
        for d in dropped:
            by.setdefault(d["reason"], []).append(d["id"])
        issues.append(_issue("fixed", "rows_dropped", f"{len(dropped)} of {n} rows dropped (" +
                             "; ".join(f"{r}: {_short(v)}" for r, v in by.items()) + ")", ids=[d["id"] for d in dropped]))
    if imputed:
        issues.append(_issue("fixed", "cells_imputed", f"{len(imputed)} cells imputed with {imputer} "
                             f"in {len({c['id'] for c in imputed})} rows"))
    rows = {"n_in": n, "n_used": int(keep.sum()), "dropped": dropped, "imputed": imputed}
    return rows, keep, Vf


def _check_args(consumer, imputer, nstds):
    if consumer not in CONSUMERS:
        raise ArgumentError(f"consumer must be one of {CONSUMERS}, not {consumer!r}")
    if imputer not in IMPUTERS:
        raise ArgumentError(f"imputer must be one of {IMPUTERS}, not {imputer!r}")
    if nstds is not None:
        nstds = float(nstds)
        if not nstds > 0:
            raise ArgumentError(f"nstds must be a positive number or None, not {nstds!r}")
    return nstds


# ---------------------------------------------------------------------------------------------------- network


def _check_network(path, role_map, first_energy, issues, declared=None):
    """Validate a rxn_network table against a role map (None in kinetic mode); `declared` is {species: role} from
    the roles file. Returns the report section."""
    declared = declared or {}
    raw, body, fixes, err = _read_table(path)
    section = {"path": str(path), "n_steps": None, "species": [], "targets": []}
    if err:
        err["message"] = "rxn_network: " + err["message"]
        issues.append(err)
        return section
    for f in fixes:
        f["message"] = "rxn_network: " + f["message"]
    issues.extend(fixes)
    hdr, ok = _headers(raw, issues, what="rxn_network")
    labels = [str(v).strip() for v in body[0].tolist()]
    low = [lab.lower() for lab in labels]
    c0_row = None
    if low and low[-1] in C0_LABELS:
        c0_row = len(low) - 1
    elif any(lab in C0_LABELS for lab in low):
        issues.append(_issue("error", "c0_not_last", "rxn_network: the c0 (initial concentrations) row must be the "
                             "last row"))
    else:
        issues.append(_issue("error", "no_c0", "rxn_network: no c0 / initial_conc row (initial concentrations) as "
                             "the last row"))
    steps = [i for i in range(len(labels)) if i != c0_row]
    section["n_steps"] = len(steps)
    if not steps:
        issues.append(_issue("error", "no_steps", "rxn_network has no elementary steps"))
    ints = {c for c, r in (role_map or {}).items() if r == "INT"}
    tss = {c for c, r in (role_map or {}).items() if r == "TS"}
    species = []
    for j in ok:
        name = hdr[j]
        conv = pd.to_numeric(body[j], errors="coerce")
        bad = body[j].notna().to_numpy() & conv.isna().to_numpy()
        vals = conv.fillna(0.0).to_numpy(dtype=float)
        if bad.any():
            issues.append(_issue("error", "non_numeric", f"rxn_network: species {name!r} has non-numeric cells "
                                 f"in rows {_short(labels[i] for i in np.where(bad)[0])}", species=name))
        coeffs = vals[steps]
        if np.any(np.abs(coeffs - np.round(coeffs)) > 1e-9):
            issues.append(_issue("error", "non_integer", f"rxn_network: species {name!r} has non-integer "
                                 "stoichiometric coefficients", species=name))
        c0 = float(vals[c0_row]) if c0_row is not None else None
        if c0 is not None and c0 < 0:
            issues.append(_issue("error", "negative_c0", f"rxn_network: species {name!r} has a negative c0",
                                 species=name))
        if role_map is not None and name in tss:
            role = None
            issues.append(_issue("error", "ts_species", f"rxn_network species {name!r} is a TS column of the "
                                 "reaction table; transition states are not species", species=name))
        elif name in declared:
            role = declared[name]
        elif role_map is not None and name in ints:
            role = "intermediate"
        elif _REACTANT_RE.match(name):
            role = "reactant"
        elif "prod" in name.lower() or _PRODUCT_RE.match(name):
            role = "target" if "*" in name else "product"
        elif role_map is None:
            role = "intermediate"
        else:
            role = None
            issues.append(_issue("error", "unmatched_species",
                                 f"rxn_network species {name!r} is not an INT column of the reaction table and has no "
                                 "reactant (R-...) or product (P-..., ...prod...) name; rename it to match the table, "
                                 "or declare its role under \"species\" in the roles file", species=name))
        if role is not None and not np.any(coeffs != 0):
            issues.append(_issue("warning", "spectator", f"rxn_network species {name!r} takes part in no step",
                                 species=name))
            role = "spectator"
        species.append({"name": name, "role": role, "c0": _num(c0) if c0 is not None else None})
    if species:
        first = species[0]
        if first["role"] == "intermediate":
            first["role"] = "catalyst"
        if role_map is not None and first_energy is not None and first["name"] != first_energy:
            issues.append(_issue("error", "catalyst_mismatch",
                                 f"rxn_network species 0 ({first['name']!r}) is the catalyst and must be the table's "
                                 f"first energy column ({first_energy!r})", species=first["name"]))
        elif role_map is None and first["role"] not in ("catalyst", "spectator"):
            issues.append(_issue("error", "catalyst_mismatch", f"rxn_network species 0 ({first['name']!r}) must be "
                                 "the catalyst, not a reactant or product", species=first["name"]))
    names_net = {s["name"] for s in species}
    for name in declared:
        if name not in names_net:
            issues.append(_issue("error", "roles_unknown_species", f"roles file declares species {name!r}, which is "
                                 "not a species of the rxn_network", species=name))
    if role_map is not None:
        for c in sorted(ints - names_net, key=list(role_map).index):
            issues.append(_issue("error", "int_not_in_network", f"INT column {c!r} of the reaction table is not a "
                                 "species of the rxn_network", column=c))
    if species and not any(s["role"] == "target" for s in species):
        issues.append(_issue("warning", "no_target", "rxn_network has no target species (a product whose name holds '*')"))
    section["species"] = species
    section["targets"] = [s["name"] for s in species if s["role"] == "target"]
    return section


# ---------------------------------------------------------------------------------------------------- profile


def _report(path, consumer, imputer, seed, nstds):
    return {"schema": SCHEMA, "clear": False, "consumer": consumer, "table": str(path), "id_column": None,
            "issues": [], "role_map": {}, "roles_source": LEGACY, "pathways": [],
            "rows": {"n_in": 0, "n_used": 0, "dropped": [], "imputed": []},
            "imputer": {"strategy": imputer, "seed": seed}, "nstds": nstds}


def _finish(rep):
    errors = [i for i in rep["issues"] if i["level"] == "error"]
    rep["clear"] = not errors and rep["rows"]["n_used"] > 0
    return rep


def load_profile(path, consumer="volcanic", imputer="none", seed=None, nstds=None, network=None, roles=None):
    """Validate a reaction table and return (report, frame). See the module docstring."""
    nstds = _check_args(consumer, imputer, nstds)
    rep = _report(path, consumer, imputer, seed, nstds)
    issues = rep["issues"]
    raw, body, fixes, err = _read_table(path)
    if err:
        issues.append(err)
        if network is not None:
            rep["network"] = _check_network(network, None, None, [])
        return _finish(rep), None
    issues.extend(fixes)
    hdr, ok = _headers(raw, issues)
    rep["id_column"] = hdr[0]
    ids = _ids(body, issues)

    declared, source, role_issues, declared_species = _load_roles(roles, path)
    rep["roles_source"] = source
    issues.extend(role_issues)
    roles_by_j = {}
    for j in ok:
        if declared is not None and hdr[j] in declared:
            roles_by_j[j] = declared[hdr[j]]
        else:
            roles_by_j[j] = role_of(hdr[j])
            if declared is not None:
                issues.append(_issue("warning", "role_by_name", f"column {hdr[j]!r} is not in the roles file; the "
                                     f"legacy name rule reads it as {roles_by_j[j]}", column=hdr[j]))
    if declared is not None:
        for name, role in declared.items():
            if name == hdr[0] and role == "ID":
                continue
            if name == hdr[0] or name not in hdr[1:]:
                issues.append(_issue("error", "roles_unknown_column", f"roles file names column {name!r}, which is "
                                     "not a data column of the table", column=name))
            elif role == "ID":
                issues.append(_issue("error", "bad_role", f"only column 0 can be the ID column, not {name!r}",
                                     column=name))
                for j in ok:
                    if hdr[j] == name:
                        roles_by_j[j] = role_of(name)
    role_map = {}
    for j in ok:
        role_map.setdefault(hdr[j], roles_by_j[j])
    rep["role_map"] = role_map
    rep["pathways"] = pathways(role_map)

    V = _numeric(body, ok, hdr, ids, issues)
    names = [hdr[j] for j in ok]
    energy_k = [k for k, j in enumerate(ok) if roles_by_j[j] in ENERGY_ROLES]
    desc_k = [k for k, j in enumerate(ok) if roles_by_j[j] == "DESCRIPTOR"]
    _structure(rep, ids, V, names, [roles_by_j[j] for j in ok], energy_k, consumer)

    rows, keep, Vf = _curate(ids, V, energy_k, desc_k, names, imputer, seed, nstds, issues)
    rep["rows"] = rows
    if rows["n_in"] and rows["n_used"] == 0:
        issues.append(_issue("error", "no_rows", f"no row left to use: all {rows['n_in']} rows were dropped"
                             + (" (imputer none: every row has a missing value; pass an imputer to fill them)"
                                if imputer == "none" else "")))
    elif not rows["n_in"]:
        issues.append(_issue("error", "no_rows", "the table has no data rows"))
    elif consumer == "volcanic" and rows["n_used"] < 3:
        issues.append(_issue("warning", "few_rows", f"only {rows['n_used']} rows left; volcanic's regressions need "
                             "several catalysts"))

    if network is not None:
        first = names[energy_k[0]] if energy_k else None
        rep["network"] = _check_network(network, role_map, first, issues, declared_species)

    frame = pd.DataFrame(Vf[keep], columns=names)
    frame.insert(0, hdr[0] or "id", [ids[i] for i in range(len(ids)) if keep[i]])
    return _finish(rep), frame


def _structure(rep, ids, V, names, roles, energy_k, consumer):
    """Reference state and the end columns each consumer needs."""
    issues = rep["issues"]
    if not energy_k:
        issues.append(_issue("error", "no_energy_columns", "the table has no energy columns (INT, TS or END)"))
        return
    k0 = energy_k[0]
    if roles[k0] != "INT":
        issues.append(_issue("error", "reference_not_int", f"the first energy column {names[k0]!r} is read as "
                             f"{roles[k0]}; it must be the reference intermediate", column=names[k0]))
    col = V[:, k0]
    bad = [i for i in range(len(col)) if not col[i] == 0.0]
    if bad:
        issues.append(_issue("error", "reference_not_zero",
                             f"the first energy column {names[k0]!r} is the reference state and must be exactly 0.0 "
                             f"in every row; {len(bad)} rows are not: "
                             f"{_short(f'{ids[i]}={_num(col[i])}' for i in bad)}",
                             column=names[k0], ids=[ids[i] for i in bad]))
    ends = [names[k] for k in energy_k if roles[k] == "END"]
    last = names[energy_k[-1]]
    if consumer == "volcanic":
        if len(ends) != 1:
            issues.append(_issue("error", "end_count",
                                 f"volcanic needs exactly one END column (the reaction energy), found {len(ends)}"
                                 + (f": {', '.join(ends)}; a multi-pathway table is for mikimo" if ends else
                                    "; name it 'Product' or declare it in a roles file"), columns=ends))
        elif ends[0] != last:
            issues.append(_issue("error", "end_not_last", f"the END column {ends[0]!r} must be the last energy column "
                                 f"(found {last!r} after it)", column=ends[0]))
    else:
        if not ends:
            issues.append(_issue("error", "end_count", "mikimo needs an END column at the end of each pathway; found "
                                 "none", columns=[]))
        elif roles[energy_k[-1]] != "END":
            issues.append(_issue("error", "end_not_last", f"the last pathway has no END: the last energy column "
                                 f"{last!r} is read as {roles[energy_k[-1]]}", column=last))


def _structure_errors(names, roles, values, ids, consumer):
    """Reference state and end-column errors for a table a consumer has already read itself (names, roles, values, ids)."""
    rep = {"issues": []}
    roles = list(roles)
    energy_k = [k for k, r in enumerate(roles) if r in ENERGY_ROLES]
    _structure(rep, [str(i) for i in ids], np.asarray(values, dtype=float), [str(n).strip() for n in names], roles,
               energy_k, consumer)
    return [i for i in rep["issues"] if i["level"] == "error"]


def validate_profile(path, consumer="volcanic", imputer="none", seed=None, nstds=None, network=None, roles=None):
    """Validate a reaction table against reaction_profile/v1. See the module docstring."""
    return load_profile(path, consumer=consumer, imputer=imputer, seed=seed, nstds=nstds, network=network,
                        roles=roles)[0]


# ---------------------------------------------------------------------------------------------------- kinetic


def validate_kinetic_data(path, network=None, imputer="none", seed=None, nstds=None, units="log10_k"):
    """Validate a mikimo kinetic_data table. See the module docstring."""
    nstds = _check_args("mikimo", imputer, nstds)
    rep = _report(path, "mikimo", imputer, seed, nstds)
    del rep["pathways"]
    rep["roles_source"] = "kinetic column names"
    issues = rep["issues"]
    kin = {"units": units, "n_steps": None, "forward": [], "reverse": [], "descriptors": [], "legacy_names": False}
    rep["kinetic"] = kin
    if units not in UNITS:
        issues.append(_issue("error", "units", f"kinetic data units must be one of {UNITS} (log10 of the rate "
                             f"constants), not {units!r}"))
    raw, body, fixes, err = _read_table(path)
    if err:
        issues.append(err)
        if network is not None:
            rep["network"] = _check_network(network, None, None, [])
        return _finish(rep)
    issues.extend(fixes)
    hdr, ok = _headers(raw, issues)
    rep["id_column"] = hdr[0]
    ids = _ids(body, issues)
    roles_by_j, fwd, rev, legacy = {}, {}, {}, []
    for j in ok:
        role, idx, leg = _kinetic_role(hdr[j])
        if role is None:
            issues.append(_issue("error", "unknown_column", f"kinetic data column {hdr[j]!r} is neither a rate "
                                 "constant (kf_<i>, kr_<i>) nor a Descriptor", column=hdr[j]))
            continue
        roles_by_j[j] = role
        if leg:
            legacy.append(hdr[j])
        if role in ("KF", "KR"):
            book = fwd if role == "KF" else rev
            if idx in book:
                issues.append(_issue("error", "duplicate_rate", f"two columns hold {role.lower()}_{idx}: "
                                     f"{book[idx]!r} and {hdr[j]!r}", columns=[book[idx], hdr[j]]))
            else:
                book[idx] = hdr[j]
    if legacy:
        kin["legacy_names"] = True
        issues.append(_issue("fixed", "legacy_rate_names", f"legacy rate-constant names read as kf_i (k_i) and kr_i "
                             f"(k_-i): {_short(legacy)}"))
    n = max(list(fwd) + list(rev) + [0])
    for book, tag in ((fwd, "kf"), (rev, "kr")):
        miss = [i for i in range(1, n + 1) if i not in book]
        if miss:
            issues.append(_issue("error", "rate_gap", f"kinetic data has no {', '.join(f'{tag}_{i}' for i in miss)} "
                                 f"(forward and reverse must both cover steps 1..{n})"))
    if n == 0:
        issues.append(_issue("error", "no_rates", "kinetic data has no rate-constant columns"))
    kin["n_steps"] = n
    kin["forward"] = [fwd[i] for i in range(1, n + 1) if i in fwd]
    kin["reverse"] = [rev[i] for i in range(1, n + 1) if i in rev]
    rate_js = [j for j in ok if roles_by_j.get(j) in ("KF", "KR")]
    if [hdr[j] for j in rate_js] != kin["forward"] + kin["reverse"]:
        issues.append(_issue("warning", "rate_order", "rate-constant columns are not in the order kf_1..kf_n, "
                             "kr_1..kr_n; read them by name, a positional reader would misassign them"))
    role_map = {hdr[j]: roles_by_j[j] for j in ok if j in roles_by_j}
    rep["role_map"] = role_map
    kin["descriptors"] = [c for c, r in role_map.items() if r == "DESCRIPTOR"]
    used = [j for j in ok if j in roles_by_j]
    V = _numeric(body, used, hdr, ids, issues)
    names = [hdr[j] for j in used]
    rate_k = [k for k, j in enumerate(used) if roles_by_j[j] in ("KF", "KR")]
    desc_k = [k for k, j in enumerate(used) if roles_by_j[j] == "DESCRIPTOR"]
    if rate_k and np.nanmax(V[:, rate_k], initial=-np.inf) > 20:
        issues.append(_issue("warning", "units_suspect", "some values exceed 20 (k > 1e20, beyond the transition-state "
                             "limit); kinetic data holds log10 k, so check that the values are not raw rate constants"))
    rows, keep, _ = _curate(ids, V, rate_k, desc_k, names, imputer, seed, nstds, issues)
    rep["rows"] = rows
    if rows["n_in"] and rows["n_used"] == 0:
        issues.append(_issue("error", "no_rows", f"no row left to use: all {rows['n_in']} rows were dropped"))
    elif not rows["n_in"]:
        issues.append(_issue("error", "no_rows", "the table has no data rows"))
    if network is not None:
        rep["network"] = _check_network(network, None, None, issues)
        ns = rep["network"]["n_steps"]
        if ns is not None and n and ns != n:
            issues.append(_issue("error", "step_count", f"kinetic data has {n} steps (kf/kr pairs) but the rxn_network "
                                 f"has {ns}"))
    return _finish(rep)


# ---------------------------------------------------------------------------------------------------- CLI


def _summary(rep):
    lines = [f"{rep['table']}: {'CLEAR' if rep['clear'] else 'REFUSED'} ({rep['consumer']}, roles: {rep['roles_source']})"]
    if rep.get("role_map"):
        lines.append("roles: " + ", ".join(f"{c}={r}" for c, r in rep["role_map"].items()))
    r = rep["rows"]
    lines.append(f"rows: {r['n_used']} of {r['n_in']} used, {len(r['dropped'])} dropped, {len(r['imputed'])} cells imputed")
    for d in r["dropped"]:
        lines.append(f"  dropped {d['id']}: {d['reason']}" + (f" (missing {', '.join(d['missing'])})" if d["missing"] else ""))
    for c in r["imputed"]:
        lines.append(f"  imputed {c['id']}/{c['column']} = {c['value']} ({c['reason']})")
    for i in rep["issues"]:
        lines.append(f"{i['level'].upper()}: {i['message']}")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m navicat_volcanic.profile",
                                 description="Validate a reaction table against the reaction_profile/v1 contract.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="validate a reaction (or kinetic) table")
    v.add_argument("table")
    v.add_argument("--consumer", default="volcanic", choices=CONSUMERS)
    v.add_argument("--network", default=None, help="mikimo rxn_network file")
    v.add_argument("--imputer", default="none", choices=IMPUTERS)
    v.add_argument("--nstds", type=float, default=None, help="outlier filter width in standard deviations (default: none)")
    v.add_argument("--seed", type=int, default=None)
    v.add_argument("--roles", default=None, help="roles file (default: <table name without its extension>.roles.json or roles.json next to "
                   "the table, else the legacy name rule)")
    v.add_argument("--kinetic", action="store_true", help="the table is mikimo kinetic_data (log10 k)")
    v.add_argument("--write-roles", default=None, help="write the resolved role map to this file")
    v.add_argument("--json", action="store_true", help="print the report as one JSON object")
    a = ap.parse_args(argv)
    try:
        if a.kinetic:
            rep = validate_kinetic_data(a.table, network=a.network, imputer=a.imputer, seed=a.seed, nstds=a.nstds)
        else:
            rep = validate_profile(a.table, consumer=a.consumer, imputer=a.imputer, seed=a.seed, nstds=a.nstds,
                                   network=a.network, roles=a.roles)
        if a.write_roles:
            save_roles(rep, a.write_roles)
    except ArgumentError as e:
        msg = f"invalid argument: {e}"
        if a.json:
            print(json.dumps({"schema": SCHEMA, "clear": False, "issues": [_issue("error", "argument", msg)]}))
        else:
            print(msg, file=sys.stderr)
        return 2
    except Exception as e:  # never a traceback on the command line
        msg = f"internal error: {type(e).__name__}: {e}"
        if a.json:
            print(json.dumps({"schema": SCHEMA, "clear": False, "issues": [_issue("error", "internal", msg)]}))
        else:
            print(msg, file=sys.stderr)
        return 3
    if a.json:
        print(json.dumps(rep, ensure_ascii=True))
    else:
        text = _summary(rep)
        try:
            print(text)
        except UnicodeEncodeError:
            print(text.encode("ascii", "backslashreplace").decode("ascii"))
    return 0 if rep["clear"] else 1


if __name__ == "__main__":
    sys.exit(main())
