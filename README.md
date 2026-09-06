# Chiral perfect-absorption design search

Code accompanying supplementary Note S5 of *Chiral Perfect Absorption in Flat Optics: A Jones Exceptional Point of Maximal Response Strength*.

## Dependencies

Required:

- [S4](https://web.stanford.edu/group/fan/S4/install.html) with its Lua frontend — the RCWA solver.
- Lua development libraries and a C/C++ compiler — needed to build S4.
- Python 3.11 or later.
- NumPy and SciPy — versions are specified in `requirements.txt`.
- Git — for cloning this repository.

Optional:

- [AAA](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.interpolate.AAA.html) (`scipy.interpolate.AAA`) — used only for refitting the TiO₂ material model; it is included in SciPy 1.15 or later, so no separate package installation is needed.

## Setup

Use an environment where the dependencies listed above are already installed.

1. Clone this repository and enter its directory.

   ```sh
   git clone https://github.com/Jhnnn001/chiralPA-design-search.git
   cd chiralPA-design-search
   ```

2. Check Python, NumPy, and SciPy in your current environment.
   Compare the reported versions with the requirements above and in `requirements.txt`.

   ```sh
   python --version
   python -c "import numpy, scipy; print('NumPy:', numpy.__version__); print('SciPy:', scipy.__version__)"
   ```

3. Run the package checks, including the S4 Lua frontend.

   ```sh
   python check.py --s4 S4
   ```

   If S4 is not on your PATH, replace `--s4 S4` with `--s4 /path/to/S4` in the check and search commands.
   Lua development libraries and a C/C++ compiler are only needed when building S4, so they are not checked here.

4. Optional: check AAA support for material refitting.

   ```sh
   python -c "from scipy.interpolate import AAA; print('AAA is available')"
   ```

## Run

Run the search from the repository directory.
The general form is:

```sh
python search.py <case> --seed <int> --trials <int> --workers <int> --numg <int> --budget-b <int> --budget-c <int>
```

The target is required and every `--` option may be omitted; defaults are given in parentheses.

- `<case>` — the target: `plain`, `nilpotent`, or `maximal`.
- `--seed <int>` — random seed; a given seed reproduces the same search exactly (default 0).
- `--trials <int>` — maximum number of independent trials; the search stops at the first accepted design (default 10).
- `--workers <int>` — parallel S4 processes for Stage A; Stages B and C evaluate one point at a time (default 4).
- `--numg <int>` — Fourier truncation order of the RCWA solver, used by every stage (default 101).
- `--budget-b <int>` — S4 evaluations allowed for each Stage B candidate (default 400).
- `--budget-c <int>` — S4 evaluations allowed for each Stage C candidate, shared by all of its starts (default 4000).
- `--population <int>` — Stage A population size, at least 3 (default 100).
- `--generations <int>` — Stage A generations (default 100).
- `--starts <int>` — Stage C starts per candidate (default 4).
- `--s4 <path>` — S4 Lua executable when it is not on your PATH (default `S4`).
- `--start <file>` — refine a saved design from Stage B onward; the file is JSON with a `p_nm` object as written in the result records.
- `--output <dir>` — output directory instead of the timestamped default.

For example:

```sh
python search.py maximal --seed 0 --trials 500 --workers 8 --numg 51
```

Each run creates one directory under `runs/`, named by the target and the UTC start time:

```
runs/<case>_<YYYYMMDD>T<HHMMSS><ffffff>Z/
├── settings.json     # the settings used, with the solver name and the NumPy and SciPy versions
├── progress.log      # UTC-stamped stage messages, appended while the search runs; follow with tail -f
├── candidates.jsonl  # one line per trial's Stage A summary and per visited Stage B and Stage C candidate
└── result.json       # found or not_found, trials and evaluations used, and the accepted design with its metrics
```

`result.json` reports `found` only after the final gate and three matching fresh evaluations; an exhausted search reports `not_found` and exits with code 1.
Example runs for the three targets are kept under `runs/`.

Other commands:

Full option list:

```sh
python search.py --help
```

Short execution check, with no expectation of finding an EP:

```sh
python search.py maximal --population 3 --generations 1 --trials 1 --budget-b 24 --budget-c 24
```

TiO₂ refit with six pole pairs and fourteen starts; it prints coefficients without modifying the frozen model used by the search:

```sh
python materials/fit_tio2.py 6 14
```

## Files and folders

| Location | Contents |
| --- | --- |
| `search.py` | Search entry point and Stages A, B, and C. |
| `geometry.py`, `objectives.py` | Geometry constraints and repair, objectives, and acceptance gates. |
| `solver.py`, `solver.lua` | S4 execution and reflection Jones-matrix extraction. |
| `materials/` | Frozen Ag, SiO₂, and TiO₂ models and the TiO₂ fitting script. |
| `materials/data/` | Johnson–Christy Ag data, Malitson SiO₂ data, and the measured TiO₂ table. |
| `runs/` | Generated settings, progress logs, candidates, and results; the three example runs are included. |
| `check.py` | Runnable numerical and execution checks. |

Material sources are identified in the model and table headers.
The Ag and SiO₂ models use published coefficients; only TiO₂ is fitted locally.
