"""Refit TiO₂ using the project AAA seed and least-squares polish.

Method follows Betz et al., Laser Photonics Rev. 18, 2400584 (2024).
Run: python materials/fit_tio2.py [pole_pairs] [starts].
Prints coefficients; the search uses the frozen coefficients in models.py.
"""
from pathlib import Path
import sys
import warnings

import numpy as np
from scipy.interpolate import AAA
from scipy.optimize import least_squares

from models import lorentz

MEASURED = Path(__file__).parent / "data" / "tio2_measured.csv"
C_NM_S = 2.99792458e17          # nm/s
W_UNIT = 1e15                   # rad/s per fitted unit

BMIN = 0.05     # damping floor; keeps AAA's near-real spurious poles out
PEN = 30.0      # weight of the Im n >= 0 penalty (measured k is >= 0 everywhere)


def load():
    d = np.genfromtxt(MEASURED, delimiter=",")
    wl, n, k = d.T
    w = 2 * np.pi * C_NM_S / wl / W_UNIT
    o = np.argsort(w)
    return wl[o], w[o], n[o], k[o]


WL, W, N, K = load()
EPS = (N + 1j * K) ** 2


def model(w, einf, Om, sig):
    """ε at (complex) ω given in units of 10¹⁵ rad/s."""
    return lorentz(w, einf, Om, sig)


def unpack(p, L):
    einf = p[0]
    a, b, sr, si = (p[1 + i * L:1 + (i + 1) * L] for i in range(4))
    return einf, a - 1j * b, sr + 1j * si


def aaa_seed(L):
    """L lower-half-plane pole pairs from AAA on the conjugate-mirrored data."""
    x = np.concatenate([-W[::-1], W])
    y = np.concatenate([np.conj(EPS[::-1]), EPS])
    for max_terms in range(2 * L + 1, 2 * L + 30, 2):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            r = AAA(x, y, rtol=1e-13, max_terms=max_terms)
        pol, res = np.array(r.poles()), np.array(r.residues())
        keep = (pol.real > 0) & (pol.imag < -BMIN) & (np.abs(res) > 1e-3)
        if keep.sum() >= L:
            idx = np.argsort(-np.abs(res[keep]))[:L]
            return float(np.mean(EPS.real)), pol[keep][idx], -1j * res[keep][idx]
    band = np.linspace(W.min(), W.max(), L + 2)[1:-1]
    return float(np.mean(EPS.real)), band - 0.5j, np.full(L, 1.0 + 0j)


def polish(L, p0):
    lo = np.concatenate([[1.0], np.full(L, 0.1), np.full(L, BMIN), np.full(2 * L, -np.inf)])
    hi = np.concatenate([[10.], np.full(L, 60.), np.full(L, 60.), np.full(2 * L, np.inf)])

    def resid(p):
        m = np.sqrt(model(W, *unpack(p, L)))
        return np.concatenate([m.real - N, m.imag - K, PEN * np.minimum(m.imag, 0.0)])

    s = least_squares(resid, np.clip(p0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi),
                      x_scale="jac", max_nfev=3000)
    einf, Om, sig = unpack(s.x, L)
    m = np.sqrt(model(W, einf, Om, sig))
    return np.abs(m.real - N).max(), np.abs(m.imag - K).max(), m.imag.min(), einf, Om, sig


def fit(L, n_starts=14, seed=0):
    """Multi-start polish; keeps the run with the smallest max|Δn|."""
    rng = np.random.default_rng(seed)
    einf0, Om0, sig0 = aaa_seed(L)
    base = np.concatenate([[einf0], Om0.real, -Om0.imag, sig0.real, sig0.imag])
    best = None
    for t in range(n_starts):
        if t == 0:
            p0 = base
        elif t <= n_starts // 2:
            p0 = base * (1 + 0.25 * rng.standard_normal(base.shape))
        else:   # cold start: poles spread over and beyond the sampled band
            p0 = np.concatenate([[2.4], rng.uniform(W.min() * 0.7, W.max() * 1.4, L),
                                 rng.uniform(0.1, 3.0, L), rng.standard_normal(2 * L)])
        try:
            r = polish(L, p0)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if best is None or r[0] < best[0]:
            best = r
    if best is None:
        raise RuntimeError("No material fit converged")
    return best


def report(L, res):
    dn, dk, kmin, einf, Om, sig = res
    m = np.sqrt(model(W, einf, Om, sig))
    print(f"L = {L} pole pairs    max|dn| = {dn:.3e}   max|dk| = {dk:.3e}   "
          f"min k = {kmin:+.2e}")
    for lo, hi in ((238, 300), (300, 400), (400, 700), (700, 1034)):
        s = (WL >= lo) & (WL < hi)
        print(f"    {lo:4d}-{hi:4d} nm : max|dn| = {np.abs(m.real - N)[s].max():.2e}   "
              f"max|dk| = {np.abs(m.imag - K)[s].max():.2e}")
    print(f"\nTIO2_EPS_INF = {einf:.8f}\nTIO2_OMEGA = np.array([")
    order = np.argsort(Om.real)
    for O in Om[order]:
        print(f"    {O.real:.8f} - {-O.imag:.8f}j,")
    print("])\nTIO2_SIGMA = np.array([")
    for S in sig[order]:
        print(f"    {S.real:+.8f} {'+' if S.imag >= 0 else '-'} {abs(S.imag):.8f}j,")
    print("])")


if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    n_starts = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    if L < 1 or n_starts < 1:
        raise ValueError("Pole-pair count and number of starts must be positive")
    report(L, fit(L, n_starts))
