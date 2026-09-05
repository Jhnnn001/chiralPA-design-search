"""S4 subprocess evaluations with geometry checks, retries, caching, and budgets."""
import os
from pathlib import Path
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from geometry import feasible, rectangles
from materials.models import permittivities
from objectives import metrics

WORKER = Path(__file__).with_suffix(".lua")


class BudgetReached(Exception):
    pass


def solve(s4, numg, p):
    eps = permittivities(p[10])
    values = np.r_[numg, p*1e-3, np.column_stack([eps.real, eps.imag]).ravel(),
                   rectangles(p).ravel()*1e-3]
    request = ";".join(format(v, ".17g") for v in values)+"\n"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    try:
        result = subprocess.run([s4, str(WORKER)], input=request, text=True,
                                capture_output=True, timeout=300, env=env)
    except subprocess.TimeoutExpired:
        return None
    if result.returncode:
        return None
    for line in reversed(result.stdout.splitlines()):
        if line.startswith("J "):
            try:
                values = np.array([float(x) for x in line.split()[1:]])
            except ValueError:
                return None
            if values.size == 8:
                return metrics((values[::2]+1j*values[1::2]).reshape(2, 2))
    return None


class Evaluator:
    def __init__(self, s4="S4", numg=101, workers=4):
        self.s4 = shutil.which(s4)
        if self.s4 is None:
            raise FileNotFoundError("Install the S4 Lua executable or pass --s4")
        self.numg = numg
        self.pool = ThreadPoolExecutor(max_workers=workers)
        self.nfev = 0
        self.nfail = 0
        self.reset()

    def reset(self, budget=None):
        self.cache = {}
        self.limit = np.inf if budget is None else self.nfev+budget

    def batch(self, points):
        points = [np.asarray(p, float) for p in points]
        pending = {p.tobytes(): p for p in points
                   if feasible(p) and p.tobytes() not in self.cache}
        for _ in range(3):
            if not pending:
                break
            if self.nfev+len(pending) > self.limit:
                raise BudgetReached
            self.nfev += len(pending)
            results = self.pool.map(lambda p: solve(self.s4, self.numg, p), pending.values())
            failed = {}
            for (key, p), m in zip(pending.items(), results):
                if m is None:
                    self.nfail += 1
                    failed[key] = p
                else:
                    self.cache[key] = m
            pending = failed
        return [self.cache.get(p.tobytes()) if feasible(p) else None for p in points]

    def one(self, p):
        return self.batch([p])[0]

    def close(self):
        self.pool.shutdown(wait=True)
