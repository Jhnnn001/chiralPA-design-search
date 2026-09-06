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
- `--start <file>` — refine a saved design from Stage B onward; the file is JSON with a `p_nm` object as written in the result records (default: none, the search starts from Stage A).
- `--output <dir>` — output directory (default `runs/<case>_<timestamp>`).

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

`result.json` reports `found` when a design is accepted and `not_found` when every trial is exhausted.
The process exits with status 0 in the first case and 1 in the second, so a shell script can test the outcome.
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

## References

The RCWA solver is S4 [1].
Stage B uses SciPy's SLSQP and Stage C uses SciPy's trust-region-reflective least squares [2–4]; neither optimizer is implemented in this repository.
The Stage A genetic algorithm is a standard implementation of tournament selection, uniform crossover, Gaussian mutation, and elitism [5–8].
The Ag model uses the Drude–Lorentz coefficients of [9] and is checked against the table of [10]; the SiO₂ model is the Sellmeier equation of [11]; the TiO₂ model is fitted to the measured table with the AAA algorithm [12] following [13].

1. V. Liu and S. Fan, "S⁴: A free electromagnetic solver for layered periodic structures," Comput. Phys. Commun. 183, 2233–2244 (2012). https://web.stanford.edu/group/fan/S4/
2. P. Virtanen et al., "SciPy 1.0: fundamental algorithms for scientific computing in Python," Nat. Methods 17, 261–272 (2020). https://scipy.org/
3. D. Kraft, "A software package for sequential quadratic programming," DFVLR-FB 88-28, DLR (1988); used through `scipy.optimize.minimize(method="SLSQP")`. https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html
4. M. A. Branch, T. F. Coleman, and Y. Li, "A subspace, interior, and conjugate gradient method for large-scale bound-constrained minimization problems," SIAM J. Sci. Comput. 21, 1–23 (1999); used through `scipy.optimize.least_squares(method="trf")`. https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
5. J. H. Holland, *Adaptation in Natural and Artificial Systems* (University of Michigan Press, 1975).
6. D. E. Goldberg, *Genetic Algorithms in Search, Optimization, and Machine Learning* (Addison-Wesley, 1989).
7. B. L. Miller and D. E. Goldberg, "Genetic algorithms, tournament selection, and the effects of noise," Complex Systems 9, 193–212 (1995).
8. G. Syswerda, "Uniform crossover in genetic algorithms," in *Proceedings of the Third International Conference on Genetic Algorithms* (Morgan Kaufmann, 1989), pp. 2–9.
9. H. S. Sehmi, W. Langbein, and E. A. Muljarov, "Optimizing the Drude–Lorentz model for material permittivity: Method, program, and examples for gold, silver, and copper," Phys. Rev. B 95, 115444 (2017).
10. P. B. Johnson and R. W. Christy, "Optical constants of the noble metals," Phys. Rev. B 6, 4370–4379 (1972).
11. I. H. Malitson, "Interspecimen comparison of the refractive index of fused silica," J. Opt. Soc. Am. 55, 1205–1209 (1965).
12. Y. Nakatsukasa, O. Sète, and L. N. Trefethen, "The AAA algorithm for rational approximation," SIAM J. Sci. Comput. 40, A1494–A1522 (2018).
13. F. Betz, M. Hammerschmidt, L. Zschiedrich, S. Burger, and F. Binkowski, "Efficient rational approximation of optical response functions with the AAA algorithm," Laser Photonics Rev. 18, 2400584 (2024).
