"""Notched TiO₂ geometry; p = [W,L,a1,b1,a2,b2,tTiO2,tSiO2,Px,Py,lam], nm."""
import numpy as np

NAMES = ["W", "L", "a1", "b1", "a2", "b2", "tTiO2", "tSiO2", "Px", "Py", "lam"]
LO = np.array([60, 60, 20, 20, 20, 20, 60, 20, 250, 250, 500], float)
HI = np.array([460, 460, 440, 440, 440, 440, 340, 500, 630, 630, 700], float)
GMIN, RHO = 20.0, 0.90


def separation(p):
    return np.hypot(max(0, p[0] - p[2] - p[4]), max(0, p[1] - p[3] - p[5]))


def margins(p):
    W, L, a1, b1, a2, b2, _, _, Px, Py, lam = p
    return np.array([W-a1-GMIN, W-a2-GMIN, L-b1-GMIN, L-b2-GMIN,
                     separation(p)-GMIN, Px-W-GMIN, Py-L-GMIN,
                     RHO*lam-Px, RHO*lam-Py])


def feasible(p):
    p = np.asarray(p, float)
    return (p.shape == (11,) and np.isfinite(p).all()
            and np.all(p >= LO-1e-9) and np.all(p <= HI+1e-9)
            and np.all(margins(p) >= -1e-9))


def repair(p):
    p = np.asarray(p, float)
    if p.shape != (11,) or not np.isfinite(p).all():
        raise ValueError("Expected eleven finite lengths in nm")
    p = np.clip(p, LO, HI)
    p[8:10] = np.minimum(p[8:10], RHO*p[10])
    p[:2] = np.clip(p[:2], LO[:2], p[8:10]-GMIN)
    p[2:6] = np.clip(p[2:6], GMIN, np.tile(p[:2], 2)-GMIN)
    if separation(p) < GMIN:
        cuts = p[2:6].copy()
        low, high = 0.0, 1.0
        for _ in range(40):
            s = (low+high)/2
            p[2:6] = GMIN+s*(cuts-GMIN)
            if separation(p) >= GMIN:
                low = s
            else:
                high = s
        p[2:6] = GMIN+low*(cuts-GMIN)
    return p


def sample(rng):
    p = rng.uniform(LO, HI)
    p[8:10] = rng.uniform(LO[8:10], np.minimum(HI[8:10], RHO*p[10]))
    p[:2] = rng.uniform(LO[:2], np.minimum(HI[:2], p[8:10]-GMIN))
    p[2:6] = rng.uniform(GMIN, np.tile(p[:2], 2)-GMIN)
    return repair(p)


def rectangles(p):
    """Disjoint rectangles [cx, cy, half_width, half_height], in nm."""
    W, L, a1, b1, a2, b2 = p[:6]
    gx, gy = W-a1-a2, L-b1-b2
    if gx >= 0:
        rows = [[-W/2+a1/2, -b1/2, a1/2, (L-b1)/2],
                [(a1-a2)/2, 0, gx/2, L/2],
                [W/2-a2/2, b2/2, a2/2, (L-b2)/2]]
    else:
        rows = [[a1/2, L/2-b1/2, (W-a1)/2, b1/2],
                [0, (b2-b1)/2, W/2, gy/2],
                [-a2/2, -L/2+b2/2, (W-a2)/2, b2/2]]
    return np.array([row for row in rows if row[2] > 0 and row[3] > 0])


def local_box(p):
    radius = 0.05*np.maximum(np.abs(p), 100)
    radius[10] = 2
    return np.maximum(p-radius, LO), np.minimum(p+radius, HI)
