"""Three-stage EP design search; see supplementary Note S5 for the method."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares, minimize

from geometry import HI, LO, NAMES, feasible, local_box, margins, repair, sample
from objectives import CASES, fitness, gate1, gate2, gate3, priority, residual, root_error
from solver import BudgetReached, Evaluator


class Stalled(Exception):
    pass


LOG = None


def log(message):
    """Print progress and append it to progress.log once the run directory exists."""
    print(message, flush=True)
    if LOG is not None:
        with LOG.open("a") as f:
            f.write(f"{datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ} {message}\n")


def stage_a(ev, rng, case, population=100, generations=100):
    ev.reset()
    pop = np.array([sample(rng) for _ in range(population)])
    fit = np.array([fitness(m, case) for m in ev.batch(pop)])
    if np.all(fit == -1):
        raise RuntimeError("No usable initial candidates; check the S4 installation and settings")
    for generation in range(generations):
        order = np.argsort(-fit, kind="stable")
        pop, fit = pop[order], fit[order]
        children = []
        while len(children) < population-2:
            parents = []
            for _ in range(2):
                indices = rng.integers(population, size=3)
                parents.append(pop[indices[np.argmax(fit[indices])]])
            child = parents[0].copy()
            if rng.random() < 0.9:
                mask = rng.random(11) < 0.5
                child[mask] = parents[1][mask]
            mask = rng.random(11) < 0.2
            child[mask] += rng.normal(size=mask.sum())*0.12*(HI-LO)[mask]
            children.append(repair(child))
        pop = np.vstack([pop[:2], children])
        fit = np.r_[fit[:2], [fitness(m, case) for m in ev.batch(children)]]
        if (generation+1) % 10 == 0:
            log(f"A generation {generation+1}: F={fit.max():.6f}")
    return pop[np.argsort(-fit, kind="stable")]


def survivors(pop, responses, case):
    """Rank by F and apply Gate 1 before diversity selection."""
    selected = []
    ranked = sorted(zip(pop, responses), key=lambda item: fitness(item[1], case), reverse=True)
    for p, m in ranked:
        if not gate1(m, case):
            continue
        if all(np.linalg.norm((p-q)/(HI-LO)) > 0.02 for q in selected):
            selected.append(p)
        if len(selected) == 6:
            break
    return selected


def stage_b(ev, p0, case, budget=400):
    ev.reset(budget)
    best = [p0.copy(), None]
    history = []

    def objective(u):
        p = LO+np.clip(u, 0, 1)*(HI-LO)
        m = ev.one(p)
        score = fitness(m, case)
        if m is not None and score > fitness(best[1], case):
            best[:] = [p, m]
        return -score

    def callback(u):
        history.append(fitness(best[1], case))
        if len(history) > 30 and history[-1]-history[-31] < 1e-6:
            raise Stalled

    u0 = np.clip((p0-LO)/(HI-LO), 0, 1)
    try:
        objective(u0)
        minimize(objective, u0, method="SLSQP", bounds=[(0, 1)]*11,
                 constraints={"type": "ineq", "fun": lambda u: margins(LO+u*(HI-LO))},
                 callback=callback, options={"maxiter": 100, "eps": 1e-3, "ftol": 1e-12})
    except (BudgetReached, Stalled):
        pass
    return best


def stage_c(ev, p0, case, rng, budget=4000, starts=4):
    ev.reset(budget)
    lo, hi = local_box(p0)
    span = hi-lo
    best, accepted = None, None

    def objective(u, weight):
        nonlocal best, accepted
        p = lo+np.clip(u, 0, 1)*span
        m = ev.one(p)
        if gate2(m, case):
            if best is None or priority(m, case) > priority(best[1], case):
                best = (p, m)
            if gate3(m, case) and (accepted is None or priority(m, case) > priority(accepted[1], case)):
                accepted = (p, m)
        return residual(m, p, case, weight, p0[10])

    try:
        for start in range(starts):
            if start and best is None:
                break
            u = ((p0 if start == 0 else best[0])-lo)/span
            if start:
                u += rng.uniform(-0.02, 0.02, 11)
            u = np.clip(u, 0, 1)
            for weight in (0.3, 0.05):
                result = least_squares(objective, u, args=(weight,), bounds=(0, 1),
                                       method="trf", jac="2-point", diff_step=3e-4,
                                       ftol=1e-15, xtol=1e-15, gtol=1e-15, max_nfev=1500)
                u = result.x
                objective(u, weight)
    except BudgetReached:
        pass
    return accepted if accepted is not None else best


def validate(ev, p, case):
    raw = []
    if not feasible(p):
        return None
    for _ in range(3):
        ev.reset()
        m = ev.one(p)
        if not gate3(m, case):
            return None
        raw.append(m["raw"])
    tolerance = 1e-12*max(1, np.abs(raw[0]).max())
    if any(np.abs(J-raw[0]).max() > tolerance for J in raw[1:]):
        return None
    return m


def record(p, m, case):
    result = dict(case=case, p_nm=dict(zip(NAMES, p.tolist())),
                  gate1=False, gate2=False, gate3=False)
    if m is not None:
        keys = ("Ap", "Am", "C", "mean", "s1", "s2", "K", "scale", "reciprocity")
        result["metrics"] = {key: float(m[key]) for key in keys}
        result["metrics"]["F"] = float(fitness(m, case))
        q = root_error(m, case)
        result["metrics"]["root_error"] = float(q) if np.isfinite(q) else None
        result.update(gate1=bool(gate1(m, case)), gate2=bool(gate2(m, case)),
                      gate3=bool(feasible(p) and gate3(m, case)),
                      J_real=m["raw"].real.tolist(), J_imag=m["raw"].imag.tolist())
    return result


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--s4", default="S4", help="S4 Lua executable")
    parser.add_argument("--numg", type=int, default=101, help="common Fourier setting for all stages")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--population", type=int, default=100)
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--budget-b", type=int, default=400)
    parser.add_argument("--budget-c", type=int, default=4000)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--start", type=Path, help="JSON with p_nm; start at Stage B")
    parser.add_argument("--output", type=Path, help="new output directory")
    args = parser.parse_args()
    positive = (args.numg, args.trials, args.workers, args.budget_b, args.budget_c, args.starts)
    if min(positive) < 1 or args.population < 3 or args.generations < 0 or args.seed < 0:
        parser.error("Counts must be positive, population >= 3, generations and seed >= 0")
    return args


def main():
    args = arguments()
    p0 = None
    if args.start:
        data = json.loads(args.start.read_text())["p_nm"]
        p0 = np.array([data[name] for name in NAMES], float)
        if not feasible(p0):
            raise ValueError("The start geometry does not satisfy the design constraints")
    global LOG
    ev = Evaluator(args.s4, args.numg, args.workers)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = args.output or Path("runs")/f"{args.case}_{stamp}"
    try:
        output.mkdir(parents=True, exist_ok=False)
        LOG = output/"progress.log"
        settings = {k: v for k, v in vars(args).items() if k not in ("s4", "start", "output")}
        settings.update(solver=Path(args.s4).name, numpy=np.__version__, scipy=scipy.__version__)
        with (output/"settings.json").open("x") as f:
            json.dump(settings, f, indent=2)
        found = None
        trials = args.trials if p0 is None else 1
        for trial in range(trials):
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, trial]))
            log(f"Trial {trial+1}/{trials}: {args.case}")
            if p0 is None:
                pop = stage_a(ev, rng, args.case, args.population, args.generations)
                responses = ev.batch(pop)
                candidates = survivors(pop, responses, args.case)
                passed = sum(bool(gate1(m, args.case)) for m in responses)
                log(f"A: gate1={passed}/{len(pop)}, selected={len(candidates)}")
                with (output/"candidates.jsonl").open("a") as f:
                    f.write(json.dumps(dict(trial=trial+1, stage="A", population=len(pop),
                                            gate1_passed=passed, selected=len(candidates)))+"\n")
            else:
                candidates = [p0]
            for p in candidates:
                p, m = stage_b(ev, p, args.case, args.budget_b)
                with (output/"candidates.jsonl").open("a") as f:
                    f.write(json.dumps(dict(trial=trial+1, stage="B", **record(p, m, args.case)))+"\n")
                log(f"B: F={fitness(m, args.case):.6f}, gate2={bool(gate2(m, args.case))}")
                if not gate2(m, args.case):
                    continue
                candidate = stage_c(ev, p, args.case, rng, args.budget_c, args.starts)
                if candidate is None:
                    continue
                p, m = candidate
                with (output/"candidates.jsonl").open("a") as f:
                    f.write(json.dumps(dict(trial=trial+1, stage="C", **record(p, m, args.case)))+"\n")
                log(f"C: root={root_error(m, args.case):.3e}, gate3={bool(gate3(m, args.case))}")
                if gate3(m, args.case):
                    checked = validate(ev, p, args.case)
                    if checked is not None:
                        found = record(p, checked, args.case)
                        break
            if found is not None:
                break
        result = dict(status="found" if found is not None else "not_found", trials_completed=trial+1,
                      evaluations=ev.nfev, failed_evaluations=ev.nfail)
        if found is not None:
            result.update(found)
        with (output/"result.json").open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
        log(f"{result['status']}: {output}, evaluations={ev.nfev}")
        return 0 if found is not None else 1
    finally:
        ev.close()


if __name__ == "__main__":
    raise SystemExit(main())
