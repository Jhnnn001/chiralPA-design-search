# Chiral perfect-absorption design search

Code accompanying supplementary Note S5 of *Chiral Perfect Absorption in Flat Optics: A Jones Exceptional Point of Maximal Response Strength*.

## Dependencies

Required:

- [S4](https://web.stanford.edu/group/fan/S4/install.html) with its Lua frontend — the RCWA solver.
- Lua development libraries and a C/C++ compiler — needed to build S4.
- Python 3.11 or later.
- NumPy and SciPy — versions are specified in `requirements.txt`.
- Git — for cloning this repository.

Optional (material refitting only):

- [AAA](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.interpolate.AAA.html) (`scipy.interpolate.AAA`) — included in SciPy 1.15 or later; no separate package installation is needed.

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

Run these commands from the repository directory.

```sh
python search.py maximal --seed 0 --trials 10 --workers 4
python search.py nilpotent --seed 0 --trials 10 --workers 4
python search.py plain --seed 0 --trials 10 --workers 4
python search.py --help
```

Each run creates a new directory under `runs/`, containing `settings.json`, `candidates.jsonl`, and `result.json`.
`candidates.jsonl` records each trial's Stage A population size, Gate 1 pass count, and selected count, followed by the visited Stage B/C candidates with `gate1`, `gate2`, and `gate3` flags.
`result.json` records `trials_completed` and reports `found` only after Gate 3 and fresh validation pass; exhausting the search returns `not_found` and exit code 1.
Settings and tolerances are in the command-line defaults and `objectives.py`.
The common solver setting is configurable with `--numg` and recorded in `settings.json`.
Search acceptance does not establish Fourier-order convergence.
Use `--output runs/my-run` to choose a new output directory; existing directories are never overwritten.
To refine a saved design from Stage B onward, use `--start design.json`, where `design.json` contains the `p_nm` object written in a candidate or result record.
This skips Stage A and Gate 1 selection; Gates 2 and 3 still apply.

The search follows Stage A → Gate 1 → Stage B → Gate 2 → Stage C → Gate 3.
Gate 1 filters the final Stage A population before selecting at most six distinct candidates in descending F order, using a normalized separation greater than 0.02.
If none pass, the next trial starts without running Stage B.

| Target | Gate 1 response tests | Gate 2 response tests |
| --- | --- | --- |
| Maximal response | C ≥ 0.85 | C ≥ 0.98 |
| Zero eigenvalue (`nilpotent`) | C within 0.10 of 0.36; A₊ within 0.10 of 1 | C within 0.03 of 0.36; A₊ within 0.02 of 1 |
| Plain | σ₁ − σ₂ within 0.10 of 0.6; Ā within 0.10 of 0.5; C > 0 | Same targets within 0.02; C > 0 |

Here C = A₊ − A₋ and Ā = (A₊ + A₋)/2.
Gates 1 and 2 additionally require `root_error` ≤ 0.30 and ≤ 0.10 respectively for every target.
The diagnostic is √η_D with S_D ≥ 0.05 for plain EPs, and max(|tr J|, √|det J|) for the zero-eigenvalue and maximal-response targets.
These two EP screening limits are provisional settings, without a search-success calibration.
Gate 3 retains the strict response tests and requires `root_error` ≤ 0.01, K ≥ 100000, S_D > 0, σ₁ ≤ 1 + 10⁻⁸, and admissible geometry, followed by three fresh matching Jones evaluations.
The former Gates 1 and 2 are now Gates 2 and 3, with the EP proximity test added to Gate 2 and the final acceptance criteria preserved.

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
