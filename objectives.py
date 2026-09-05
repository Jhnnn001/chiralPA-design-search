"""Objectives and gates of supplementary Note S5."""
import numpy as np

CASES = ("plain", "nilpotent", "maximal")
PLUS = np.array([1, 1j])/np.sqrt(2)
MINUS = PLUS.conj()
S_MIN, Q_MAX, K_MIN = 0.05, 1e-2, 1e5
C_MIN = 0.98
# Root screening happens once, at Gate 1's level; Gate 2 tiers on the response and
# keeps the same root limit only as a guard against Stage B making the root worse.
Q_GATE1, Q_GATE2 = 0.30, 0.30
# Objective scale of each error component, kept independent of the gate limits.
Q_SCALE = 0.10
TOL = {"plain": np.array([Q_SCALE, 0.02, 0.02, 0.02]),
       "nilpotent": np.array([0.02, 0.05, Q_SCALE/2]),
       "maximal": np.array([0.02, 0.02, Q_SCALE/2])}


def metrics(raw):
    """Validate raw reciprocity before using the symmetric Jones matrix."""
    raw = np.asarray(raw, complex)
    if raw.shape != (2, 2) or not np.isfinite(raw).all():
        return None
    reciprocity = abs(raw[0, 1]-raw[1, 0])
    if reciprocity > 1e-8*max(np.linalg.norm(raw, 2), 1e-12):
        return None
    J = (raw+raw.T)/2
    a, b, d = J[0, 0], J[0, 1], J[1, 1]
    tr, det = a+d, a*d-b*b
    D = (a-d)**2+4*b*b
    scale = abs(a-d)**2+4*abs(b)**2
    s1, s2 = np.linalg.svd(J, compute_uv=False)
    if scale == 0:
        K = 1.0
    elif D == 0:
        K = 1e16
    else:
        V = np.linalg.eig(J)[1]
        V /= np.linalg.norm(V, axis=0)
        K = min(1e16, 1/max(abs(np.linalg.det(V))**2, 1e-300))
    Ap, Am = [1-np.linalg.norm(J @ v)**2 for v in (PLUS, MINUS)]
    return dict(raw=raw, J=J, tr=tr, det=det, D=D, scale=scale,
                eta=abs(D)/scale if scale else np.inf, K=K,
                s1=s1, s2=s2, Ap=Ap, Am=Am, C=Ap-Am,
                mean=(Ap+Am)/2, reciprocity=reciprocity)


def profile(m, case):
    if case == "plain":
        return np.array([m["s1"]-m["s2"]-0.6, m["mean"]-0.5])
    return np.array([m["Ap"]-1, m["Am"]-(0.64 if case == "nilpotent" else 0)])


def fitness(m, case):
    """Saturating and gate-scaled, so selection ranks by distance to the Gate 2 boundary."""
    if m is None or (case == "plain" and m["scale"] < S_MIN):
        return -1.0
    if case == "plain":
        error = np.r_[np.sqrt(m["eta"]), profile(m, case), max(0, -m["C"])]
    else:
        error = np.r_[profile(m, case), abs(m["tr"])/2]
    return 1/(1+np.linalg.norm(error/TOL[case]))


def root_error(m, case):
    if case == "plain":
        return np.sqrt(m["eta"]) if m["scale"] >= S_MIN else np.inf
    return max(abs(m["tr"]), np.sqrt(abs(m["det"])))


def residual(m, p, case, weight, lam0):
    if m is None:
        return np.full(6 if case == "plain" else 7, 6.0)
    if case == "plain":
        return np.r_[m["D"].real, m["D"].imag, weight*profile(m, case),
                     max(0, S_MIN-m["scale"]), max(0, -m["C"])]
    return np.r_[m["tr"].real, m["tr"].imag, m["det"].real, m["det"].imag,
                 weight*profile(m, case), weight*(p[10]-lam0)/lam0]


def response_gate(m, case, loose=False):
    if m is None:
        return False
    if case == "maximal":
        return m["C"] >= (0.85 if loose else C_MIN)
    tolerance = 0.10 if loose else 0.02
    if case == "nilpotent":
        return (abs(m["C"]-0.36) <= (0.10 if loose else 0.03)
                and abs(m["Ap"]-1) <= tolerance)
    return np.all(np.abs(profile(m, case)) <= tolerance) and m["C"] > 0


def gate1(m, case):
    return response_gate(m, case, loose=True) and root_error(m, case) <= Q_GATE1


def gate2(m, case):
    return response_gate(m, case) and root_error(m, case) <= Q_GATE2


def gate3(m, case):
    return (gate2(m, case) and m["scale"] > 0 and m["s1"] <= 1+1e-8
            and root_error(m, case) <= Q_MAX and m["K"] >= K_MIN)


def priority(m, case):
    error = np.abs(profile(m, case)).sum()
    if case == "plain":
        error += max(0, -m["C"])
    return -np.log10(max(root_error(m, case), 1e-16))-10*error
