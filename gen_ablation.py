#!/usr/bin/env python3
"""Generate ablation (blocking ON vs OFF) + convergence data for Figure 2."""

import numpy as np

from mc_engine import pack
from pore_accessibility import PoreGrid

FRAMEWORK = "frameworks/zeolites/FAU.cif"
STERIC = "frameworks/zeolites/FAU.blocked.npz"
GUEST = {"name": "CO2", "file": "molecules/CO2.xyz"}

LOADINGS = [8, 16, 24, 32, 40, 48]
SEEDS = [0, 1, 2, 3, 4]
MAX_ITERS = 8000

CONV_LOADING = 70          # dense loading so annealing trajectory is non-trivial
CONV_SEEDS = [0, 1, 2, 3, 4]
CONV_MAX_ITERS = 12000
CONV_RECORD = 10

OUT_ABL = "runs/ablation.npz"
OUT_CONV = "runs/convergence.npz"


def main():
    grid = PoreGrid.load(STERIC)

    # ---- ablation: misplacement vs loading, ON vs OFF ----
    mis_on = np.zeros((len(LOADINGS), len(SEEDS)))
    mis_off = np.zeros((len(LOADINGS), len(SEEDS)))
    for li, n in enumerate(LOADINGS):
        for si, seed in enumerate(SEEDS):
            r_on = pack(FRAMEWORK, [{**GUEST, "count": n}], use_blocking=True,
                        steric_grid=grid, seed=seed, max_iters=MAX_ITERS)
            r_off = pack(FRAMEWORK, [{**GUEST, "count": n}], use_blocking=False,
                         steric_grid=grid, seed=seed, max_iters=MAX_ITERS)
            mis_on[li, si] = r_on.n_misplaced / n
            mis_off[li, si] = r_off.n_misplaced / n
        print(f"[ablation] N={n}: ON={mis_on[li].mean():.3f} OFF={mis_off[li].mean():.3f}",
              flush=True)
    np.savez_compressed(OUT_ABL, loadings=np.array(LOADINGS),
                        mis_on=mis_on, mis_off=mis_off, seeds=np.array(SEEDS))
    print(f"[done] {OUT_ABL}")

    # ---- convergence traces (blocking ON, higher loading) ----
    traces_e, traces_t, traces_acc = [], [], []
    for seed in CONV_SEEDS:
        r = pack(FRAMEWORK, [{**GUEST, "count": CONV_LOADING}], use_blocking=True,
                 steric_grid=grid, seed=seed, max_iters=CONV_MAX_ITERS,
                 record_every=CONV_RECORD)
        traces_e.append(r.energy_trace)
        traces_t.append(r.temperature_trace)
        traces_acc.append(r.acceptance_trace)
        print(f"[conv] seed={seed} converged_iter={r.converged_iter} "
              f"n_rec={len(r.energy_trace)}", flush=True)
    # pad to equal length
    L = max(len(t) for t in traces_e)
    def padarr(lst):
        a = np.full((len(lst), L), np.nan)
        for i, t in enumerate(lst):
            a[i, :len(t)] = t
        return a
    np.savez_compressed(
        OUT_CONV,
        steps=np.arange(1, L + 1) * CONV_RECORD,
        energy=padarr(traces_e),
        temperature=padarr(traces_t),
        acceptance=padarr(traces_acc),
    )
    print(f"[done] {OUT_CONV}")


def conv_only():
    grid = PoreGrid.load(STERIC)
    traces_e, traces_t, traces_acc = [], [], []
    for seed in CONV_SEEDS:
        r = pack(FRAMEWORK, [{**GUEST, "count": CONV_LOADING}], use_blocking=True,
                 steric_grid=grid, seed=seed, max_iters=CONV_MAX_ITERS,
                 record_every=CONV_RECORD)
        traces_e.append(r.energy_trace)
        traces_t.append(r.temperature_trace)
        traces_acc.append(r.acceptance_trace)
        print(f"[conv] seed={seed} converged_iter={r.converged_iter} "
              f"n_rec={len(r.energy_trace)}", flush=True)
    L = max(len(t) for t in traces_e)
    def padarr(lst):
        a = np.full((len(lst), L), np.nan)
        for i, t in enumerate(lst):
            a[i, :len(t)] = t
        return a
    np.savez_compressed(
        OUT_CONV,
        steps=np.arange(1, L + 1) * CONV_RECORD,
        energy=padarr(traces_e),
        temperature=padarr(traces_t),
        acceptance=padarr(traces_acc),
    )
    print(f"[done] {OUT_CONV}")


if __name__ == "__main__":
    import sys
    if "--conv-only" in sys.argv:
        conv_only()
    else:
        main()
