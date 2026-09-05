"""Run the release checks; optionally add --s4 S4 for actual solver checks."""
import argparse
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

import search
import solver
from geometry import HI, LO, feasible, local_box, rectangles, repair, sample
from materials.models import (C_NM_S, HBAR_EVS, ag_epsilon, lorentz,
                              permittivities, sio2_epsilon, tio2_epsilon)
from objectives import (MINUS, PLUS, fitness, gate1, gate2, metrics,
                        residual, root_error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s4", help="also check the S4 Lua executable")
    args = parser.parse_args()
    rng = np.random.default_rng(3)
    for _ in range(100):
        for p in (sample(rng), repair(rng.uniform(LO-100, HI+100))):
            assert feasible(p) and np.allclose(repair(p), p)
            area = 4*np.prod(rectangles(p)[:, 2:], axis=1).sum()
            assert np.isclose(area, p[0]*p[1]-p[2]*p[3]-p[4]*p[5])
            lo, hi = local_box(p)
            assert np.all(lo <= p) and np.all(hi >= p)
            assert np.all(lo >= LO) and np.all(hi <= HI)
    for a in (100, 140):
        p = np.array([200, 200, a, 40, a, 50, 100, 100, 300, 300, 600], float)
        assert feasible(p)
        assert np.isclose(4*np.prod(rectangles(p)[:, 2:], axis=1).sum(), 40000-a*90)
    assert len(rectangles(np.array([200, 200, 100, 40, 100, 50]))) == 2
    assert not feasible(np.full(11, np.nan)) and not feasible(np.zeros(10))

    ideal = {}
    for case, alpha, beta in (("plain", np.sqrt(0.32), 0.6),
                              ("nilpotent", 0, 0.6), ("maximal", 0, 1)):
        J = alpha*np.eye(2)+beta*np.outer(PLUS, MINUS.conj())
        m = metrics(J)
        ideal[case] = m
        assert gate1(m, case) and gate2(m, case)
        assert abs(m["C"]-beta**2) < 1e-12
        assert abs(m["s1"]*m["s2"]-alpha**2) < 1e-12
        assert abs(fitness(m, case)-1) < 1e-7 and root_error(m, case) < 1e-7
        assert not gate1(metrics(J.conj()), case)
        assert np.linalg.norm(residual(m, p, case, 0.3, p[10])) < 1e-12
    scalar = metrics(0.3*np.eye(2))
    assert scalar["K"] == 1 and fitness(scalar, "plain") == -1
    assert not gate2(metrics(np.zeros((2, 2))), "nilpotent")
    assert not gate2(metrics(1.1*ideal["maximal"]["J"]), "maximal")
    assert metrics(np.array([[0, 0.1], [0.2, 0]])) is None
    perturbed = metrics(ideal["maximal"]["J"]+0.03*np.eye(2))
    assert gate1(perturbed, "maximal") and not gate2(perturbed, "maximal")
    for case, n_eq in (("plain", 2), ("nilpotent", 4), ("maximal", 4)):
        a = residual(perturbed, p, case, 0.3, p[10]-1)
        b = residual(perturbed, p, case, 0.05, p[10]-1)
        assert np.allclose(a[:n_eq], b[:n_eq])
        stop = 4 if case == "plain" else 7
        assert np.allclose(a[n_eq:stop]/6, b[n_eq:stop])
        if case == "plain":
            assert np.allclose(a[4:], b[4:])

    data = Path(__file__).parent/"materials"/"data"
    wl, n, k = np.loadtxt(data/"tio2_measured.csv", delimiter=",").T
    index = np.sqrt(tio2_epsilon(2*np.pi*C_NM_S/wl))
    assert np.max(abs(index.real-n)) < 5e-3 and np.max(abs(index.imag-k)) < 1.5e-3
    ag = np.genfromtxt(data/"ag_johnson_christy_1972_table1.csv", delimiter=",")
    ag = ag[np.isfinite(ag).all(axis=1)]
    visible = ag[(ag[:, 0] >= 1.59) & (ag[:, 0] <= 3.27)]
    index = np.sqrt(ag_epsilon(visible[:, 0]/HBAR_EVS))
    assert max(abs(index.real-visible[:, 1])) < 0.05
    assert max(abs(index.imag-visible[:, 2])) < 0.03
    si = np.genfromtxt(data/"sio2_malitson_1965_table1.csv", delimiter=",", usecols=(0, 2))
    si = si[np.isfinite(si).all(axis=1)]
    index = np.sqrt(sio2_epsilon(2*np.pi*C_NM_S/(si[:, 0]*1000)))
    assert max(abs(index.real-si[:, 1])) < 1.5e-6 and np.all(index.imag == 0)
    for wavelength in np.linspace(500, 700, 101):
        assert np.all(permittivities(wavelength).imag >= 0)
    z = np.array([3e15-0.4e15j, 4.7e15-1.1e15j])
    for model in (tio2_epsilon, ag_epsilon, sio2_epsilon):
        assert np.allclose(model(z).conj(), model(-z.conj()))
    assert np.isfinite(lorentz([2, 3], 1, [4-1j], [2+1j])).all()

    ev = solver.Evaluator(sys.executable, workers=2)
    try:
        with patch("solver.solve", return_value=ideal["maximal"]):
            ev.reset(1)
            m = ev.batch([p, p, np.zeros(11)])
            assert m[0] is m[1] and m[2] is None and ev.nfev == 1
        with patch("solver.solve", return_value=None):
            before = ev.nfev
            ev.reset(2)
            try:
                ev.one(p)
                raise AssertionError("Budget must include retries")
            except solver.BudgetReached:
                assert ev.nfev-before == 2
        for case, m in ideal.items():
            with patch("solver.solve", return_value=m):
                pop = search.stage_a(ev, rng, case, population=3, generations=1)
                p1, m1 = search.stage_b(ev, pop[0], case, budget=30)
                assert feasible(p1) and gate1(m1, case)
                with patch("search.least_squares", wraps=search.least_squares) as calls:
                    p2, m2 = search.stage_c(ev, p1, case, rng, budget=100, starts=1)
                    assert [call.kwargs["args"][0] for call in calls.call_args_list] == [0.3, 0.05]
                assert feasible(p2) and gate2(m2, case)
                before = ev.nfev
                assert search.validate(ev, p2, case) is not None
                assert ev.nfev-before == 3
    finally:
        ev.close()
    print("Geometry, materials, objectives, gates, stages, retries, and budgets: passed")

    if args.s4:
        ev = solver.Evaluator(args.s4, workers=1)
        try:
            flat = np.array([100, 150, 20, 20, 20, 20, 0, 0, 300, 300, 600], float)
            m = solver.solve(ev.s4, ev.numg, flat)
            n_ag = np.sqrt(ag_epsilon(2*np.pi*C_NM_S/flat[10]))
            assert m is not None
            assert np.allclose(m["J"], (1-n_ag)/(1+n_ag)*np.eye(2), atol=1e-10)
            m = ev.one(p)
            assert m is not None and m["s1"] <= 1+1e-8
        finally:
            ev.close()
        print("S4 Fresnel anchor and patterned Jones evaluation: passed")


if __name__ == "__main__":
    main()
