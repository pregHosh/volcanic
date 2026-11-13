volcanic: Automated Generator of Volcano Plots
==============================================
[![DOI](https://zenodo.org/badge/381737392.svg)](https://zenodo.org/badge/latestdoi/381737392)

![volcanic logo](./images/volcanic_logo.png)
[![PyPI version](https://badge.fury.io/py/navicat_volcanic.svg)](https://badge.fury.io/py/navicat_volcanic)

## Contents
* [About](#about-)
* [Install](#install-)
* [Examples](#examples-)
* [Citation](#citation-)

## About [↑](#about)

The code runs on pure python with the following dependencies: 
- `numpy`
- `scipy`
- `matplotlib`
- `pandas`
- `scikit-learn`
- `openpyxl` (to read .xlsx tables)
- Optionally, `h5py`

## Install [↑](#install)

You can install volcanic using pip:

```python
pip install navicat_volcanic
```

Afterwards, you can call volcanic as:

```python
python -m navicat_volcanic [-h] [-version] -i [FILENAMES] [-df DFILENAMES] [-nd ND] [-v VERB] [-r RUNMODE] [-lsfer | -thermo | -kinetic | -es | -tof | -all] [-T TEMP] [-pm PLOTMODE] [-ic IC] [-fc FC]
                [-rm RMARGIN] [-lm LMARGIN] [-np NPOINTS] [-d] [-is IMPUTER_STRAT] [-refill]
                [--descriptor auto|COLUMN] [--outdir OUTDIR] [--seed SEED] [--roles ROLES]
```
or simply
```python
navicat_volcanic [-h] [-version] -i [FILENAMES] [-df DFILENAMES] [-nd ND] [-v VERB] [-r RUNMODE] [-lsfer | -thermo | -kinetic | -es | -tof | -all] [-T TEMP] [-pm PLOTMODE] [-ic IC] [-fc FC]
                [-rm RMARGIN] [-lm LMARGIN] [-np NPOINTS] [-d] [-is IMPUTER_STRAT] [-refill]
                [--descriptor auto|COLUMN] [--outdir OUTDIR] [--seed SEED] [--roles ROLES]
```

Alternatively, you can download the package and execute:

```python 
python setup.py install
```

Afterwards, you can call volcanic as:

```python
python -m navicat_volcanic [-h] [-version] -i [FILENAMES] [-df DFILENAMES] [-nd ND] [-v VERB] [-r RUNMODE] [-lsfer | -thermo | -kinetic | -es | -tof | -all] [-T TEMP] [-pm PLOTMODE] [-ic IC] [-fc FC]
                [-rm RMARGIN] [-lm LMARGIN] [-np NPOINTS] [-d] [-is IMPUTER_STRAT] [-refill]
                [--descriptor auto|COLUMN] [--outdir OUTDIR] [--seed SEED] [--roles ROLES]
```
or
```python
navicat_volcanic [-h] [-version] -i [FILENAMES] [-df DFILENAMES] [-nd ND] [-v VERB] [-r RUNMODE] [-lsfer | -thermo | -kinetic | -es | -tof | -all] [-T TEMP] [-pm PLOTMODE] [-ic IC] [-fc FC]
                [-rm RMARGIN] [-lm LMARGIN] [-np NPOINTS] [-d] [-is IMPUTER_STRAT] [-refill]
                [--descriptor auto|COLUMN] [--outdir OUTDIR] [--seed SEED] [--roles ROLES]
```


Options can be consulted using the `-h` flag in either case. The help menu is quite detailed. 

To use volcanic from your own Python code, call `navicat_volcanic.volcanic.run_volcanic([...])` with the command-line arguments as a list; it returns the exit status and writes the same files as the command line. The functions that build the volcano plots and activity maps inside `volcanic.py` (`volcanic_2d`, `volcanic_3d`) are nested in the run and cannot be imported; the plotting functions they call can, from `navicat_volcanic.plotting2d` (`plot_2d_lsfer`, `plot_2d_t_volcano`, `plot_2d_k_volcano`, `plot_2d_es_volcano`, `plot_2d_tof_volcano`) and `navicat_volcanic.plotting3d` (`plot_3d_lsfer`, `plot_3d_t_volcano`, `plot_3d_k_volcano`, `plot_3d_es_volcano`, `plot_3d_tof_volcano`).

## Examples [↑](#examples)

The examples subdirectory contains a copious amount of tests which double as examples. Any of the data files can be run as:

```python
python -m navicat_volcanic -i [FILENAME]
```

At a terminal, this will query the user for options and generate the volcano plots as png images. Options can be consulted with the `-h` flag. To run without prompts, see [Running without prompts](#running-without-prompts-scripts-pipelines-agents).

The input of volcanic is a `pandas` compatible dataframe, which includes plain .csv and .xls files. 

Regarding format, volcanic expects headers for all columns. The first column must contain names/identifiers. Then, volcanic expects a number of columns with relative free energies for the species in the catalytic cycle (in order of appearance), whose headers must contain "TS" if the species is a transition state, and a final column whose header is "Product" containing the reaction energy. Non-energy descriptors can be input as a separate file using the `-df` flag or as extra columns whose headers contain the word "Descriptor".

### Checking a table before a run

Column roles are read with one name rule shared with [mikimo](https://github.com/lcmd-epfl/mikimo) (`navicat_volcanic.profile.role_of`): a header containing "descriptor" (any case) is a descriptor; "TS" in any case, not preceded by a letter and not followed by a lowercase letter, is a transition state (`TS1`, `ts2`, `ΔGRRS(TSOA)`; not `OTs` or `adducts`); a header containing "prod" (any case), starting with `P-` or `P_`, or just `P`, ends the pathway (`Product`, `P-A`; not `P2` or `Pd0L2`); anything else is an intermediate. volcanic needs exactly one end column, as the last energy column, and the first energy column (the reference state) must be exactly 0.0 in every row; tables that break these rules are refused with a message instead of being guessed at. Rows dropped for missing values or as outliers, and cells filled by an imputer, are always printed.

To check a table without running volcanic:

```
python -m navicat_volcanic.profile validate examples/example.xlsx --consumer volcanic --json
```

This prints the resolved role of every column, every refusal reason, every change applied (trimmed headers, a stripped byte-order mark), and every dropped row or imputed cell; it exits 0 when the table is usable and 1 when it is not. Roles can be declared instead of guessed in a roles file next to the table, named like the table without its extension plus `.roles.json` (`suzuki.roles.json` for `suzuki.xlsx`), or `roles.json`: `{"column header": "INT" | "TS" | "END" | "DESCRIPTOR"}`, optionally wrapped as `{"roles": {...}, "species": {"L": "reactant"}}` to also declare the role of a mikimo network species by name (`reactant`, `product`, `target`, `intermediate` or `spectator`). Missing values are not imputed unless asked (`--imputer knn|simple|iterative`), and the outlier filter runs only with `--nstds`. The full contract is in the docstring of `navicat_volcanic/profile.py`.

### Running without prompts (scripts, pipelines, agents)

```
python -m navicat_volcanic -i examples/suzuki.xlsx -es --descriptor auto --outdir results/suzuki
```

- `--descriptor auto` takes the candidate with the best mean r² of the linear scaling relations (transition-state columns are never candidates). A tie (the same r² to 4 decimals) is broken with `--seed` (default 0) and reported in the result file; a tie with a transition-state column is reported too. `--descriptor COLUMN` names the column by its header (trimmed; a `Descriptor X` column can be named `X`). With `-nd 2`, pass `--descriptor auto` or `--descriptor` twice.
- Give the volcano plots to build with a run mode: `-lsfer`, `-thermo`, `-kinetic`, `-es`, `-tof` or `-all` (or `-r 0..5`).
- When stdin is not a terminal (a script, a pipe, a scheduler), volcanic never prompts. Without `--descriptor` it stops (exit status 2) and lists the candidate descriptors with their r² and MAE, on screen and in `volcanic_result.json`; without a run mode it stops (exit status 2) naming the run-mode flags. Answers piped into stdin are not read. At a terminal, volcanic still asks as before when these options are left out.
- `--outdir DIR` puts every file the run writes in `DIR` (created when missing), including the `-d` hdf5 dump; without it files go to the working directory.
- `--seed N` also seeds the iterative imputer and volcanic's other random tie-breaks.
- `--roles FILE` declares column roles (see `roles.json` above) instead of the sidecar file or the name rule.
- `-df FILE` adds descriptor columns from a separate table: one column per descriptor, rows in the order of the reaction table (or a first column of ids, matched to the table's ids), no other columns. Each one becomes a `Descriptor <header>` column after the energy columns.

In Python, `navicat_volcanic.helpers.choose_1_dv` and `choose_2_dv` choose descriptors the same way without asking (`choice="auto"`, a column name or an index); `user_choose_1_dv` and `user_choose_2_dv` take a `choice` argument too and never call `input()` when stdin is not a terminal.

### Outputs

Every run writes `volcanic_result.json` (schema `volcanic_result/v1`), also when it stops on a refused table or a missing choice. It holds:

- `status` (`ok`, `refused`, `error`, or `declined` when no descriptor was taken at the prompt) and `error`;
- `inputs`: the input files with their sha256, `T_K`, `nd`, `runmode`, `imputer`, `nstds` (the outlier filter width, 3), `seed`;
- `validation` (the reaction_profile/v1 check of each table) and `role_map` (the role of every column);
- `rows`: `n_in`, `n_used`, every dropped row with its reason (`incomplete`, `outlier`, `not imputable`) and every imputed cell;
- `descriptor`: `chosen`, all `candidates` with `r2` and `mae`, and `how` it was chosen (ties included);
- `lsr`: for each regressed column, the `slope`, `intercept` and `r2` of its fit on the descriptor;
- `volcanoes`: per volcano kind (`thermodynamic`, `kinetic`, `energy_span`, `tof`) the `csv` and `png` file names (relative to the output directory), labels, `units` and, for 1 descriptor, the region boundaries of the curve;
- `points`: per volcano kind, one entry per catalyst `{id, x, y, determining}` (`determining` is the pair of states that sets the value: TDI and TDTS for energy span and TOF);
- `files` (everything written) and `warnings`.

The file is strict JSON: an undefined number (NaN, for example the r² of a constant column, or a point whose value cannot be computed) is written as `null`.

Exit status: 0 the run finished; 1 the input was refused (a table that breaks the contract, a bad descriptor file, no row left); 2 a choice is missing or invalid (no `--descriptor` or no run mode while stdin is not a terminal, an unknown column); 3 an internal error (the traceback is in the result file).

CSV files are written at every verbosity level and can be used to plot your volcano plot or activity map using external tools. An example is found in the `pretty_plotting_example` directory in this repository. File names are built from the column headers, with spaces, brackets, quotes and path or shell characters replaced by `_`. Higher verbosity levels (`-v 1`, `-v 2`, etc.) print more detail as the run goes. This can be useful for understanding what the code is doing and what can have possibly gone wrong. To be as automated as possible, reasonable default values are set for most choices. The generated csvs also contain the 95% confidence interval value for the plot in question, which are obtained by propagating the uncertainty of the LSFERs involved in the computation. To simplify the propagation of uncertainties in the TOF computation, the uncertainty is approximated using the energy span in TOF volcano plots. 

The plotmode (`-pm 1`, `-pm 2`) option can be used to modify the default look of the generated pngs, including more detail as the plotmode level increases. 


## Citation [↑](#citation)

Please cite the accompanying manuscript, which clarifies the details of volcano plot construction. You can find it [here](https://rdcu.be/cT7uu) and in the reference:

```
Laplaza, R., Das, S., Wodrich, M.D. et al. Constructing and interpreting volcano plots and activity maps to navigate homogeneous catalyst landscapes. Nat Protoc (2022). https://doi.org/10.1038/s41596-022-00726-2
```

Kinetic volcano plots were introduced [here](https://doi.org/10.1039/C6SC01660J) and turnover frequency/energy span volcano plots were introduced [here](https://doi.org/10.1021/acscatal.9b00717). The energy span model was introduced by Kozuch and Shaik [here](https://doi.org/10.1021/ar1000956). Please include those citations where relevant. A comprehensive account can be found [here](https://doi.org/10.1021/acs.accounts.0c00857).


---


