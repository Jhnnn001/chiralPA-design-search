# Chiral perfect-absorption design search

Code accompanying supplementary Note S5 of *Chiral Perfect Absorption in Flat Optics: A Jones Exceptional Point of Maximal Response Strength*.

## Install

Use Python 3.11 or later and the [S4 Lua executable](https://web.stanford.edu/group/fan/S4/install.html).
Install S4 separately and make the executable available as `S4` on your PATH, or supply it with `--s4`.
The upstream build requires a C/C++ compiler and Lua development libraries; its installation page lists the optional numerical libraries.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python check.py
python check.py --s4 S4
```

NumPy and SciPy are the Python dependencies; SciPy 1.15 or later supplies [AAA](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.interpolate.AAA.html) for the optional material refit.

## Run

Run these commands from the repository directory.

```sh
python search.py maximal --seed 0 --trials 10 --workers 4
python search.py nilpotent --seed 0 --trials 10 --workers 4
python search.py plain --seed 0 --trials 10 --workers 4
python search.py --help
```

Each run creates a new directory under `runs/`, containing `settings.json`, the visited Stage B/C candidates in `candidates.jsonl`, and `result.json`.
`result.json` reports `found` only after both gates and fresh validation pass; exhausting the search returns `not_found` and exit code 1.
Settings and tolerances are in the command-line defaults and `objectives.py`.
The common solver setting is configurable with `--numg` and recorded in `settings.json`.
Search acceptance does not establish Fourier-order convergence.
Use `--output runs/my-run` to choose a new output directory; existing directories are never overwritten.
To refine a saved design from Stage B onward, use `--start design.json`, where `design.json` contains the `p_nm` object written in a candidate or result record.

For a short execution check, with no expectation of finding an EP:

```sh
python search.py maximal --population 3 --generations 1 --trials 1 --budget-b 24 --budget-c 24
```

To rerun the TiO₂ fitting procedure with six pole pairs and fourteen starts:

```sh
python materials/fit_tio2.py 6 14
```

The refit prints coefficients without modifying the frozen model used by the search.

## Files and folders

| Location | Contents |
| --- | --- |
| `search.py` | Search entry point and Stages A, B, and C. |
| `geometry.py`, `objectives.py` | Geometry constraints and repair, objectives, and acceptance gates. |
| `solver.py`, `solver.lua` | S4 execution and reflection Jones-matrix extraction. |
| `materials/` | Frozen Ag, SiO₂, and TiO₂ models and the TiO₂ fitting script. |
| `materials/data/` | Johnson–Christy Ag data, Malitson SiO₂ data, and the measured TiO₂ table. |
| `runs/` | Generated settings, candidates, and results; excluded from Git. |
| `check.py` | Runnable numerical and execution checks. |

Material sources are identified in the model and table headers.
The Ag and SiO₂ models use published coefficients; only TiO₂ is fitted locally.
