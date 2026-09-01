#!/usr/bin/env python3
"""Reviewer-requested convergence, robustness, and failure diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from ovito.io import import_file

from error_metrics import evaluate, split_molecules
from gen_generalization import FRAMEWORKS, GUESTS, load_or_build_steric, target_loading
from mc_engine import pack
from pore_accessibility import PoreGrid, load_framework_ovito

OUT = Path("runs/reviewer")
OUT.mkdir(parents=True, exist_ok=True)


def write_csv(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(name):
    path = OUT / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def grid_sensitivity():
    rows = []
    for fw_name, (fw_path, family) in FRAMEWORKS.items():
        pos, cell, inv, _ = load_framework_ovito(fw_path)
        probe = 2.5 if fw_name == "FAU" else 1.8
        for n in (32, 48, 64, 96):
            for nmin in ((1, 4, 8, 16) if n == 64 else (4,)):
                g = PoreGrid.build(pos, cell, inv, grid_n=n, probe_radius=probe,
                                   min_cluster_size=nmin, framework_file=fw_path,
                                   mode="reviewer_sensitivity")
                s = g.stats
                rows.append({
                    "framework": fw_name, "family": family, "probe_A": probe,
                    "grid_n": n, "n_min": nmin,
                    "accessible_cell_fraction": s["accessible_voxels"] / n**3,
                    "accessible_void_fraction": s["accessible_fraction_of_void"],
                    "blocked_void_fraction": s["blocked_voxels"] / max(1, s["void_voxels"]),
                    "blocked_clusters": s.get("n_clusters", 0),
                    "largest_blocked_cluster": (s.get("largest_clusters") or [0])[0],
                    "speckle_removed": s.get("speckle_removed", 0),
                })
                print("grid", fw_name, n, nmin, flush=True)
    write_csv("grid_and_filter_sensitivity.csv", rows)
    return rows


def run_case(fw_name, guest_name, cutoff, penalty, seed, max_iters):
    fw_path, _ = FRAMEWORKS[fw_name]
    guest_file, guest_radius = GUESTS[guest_name]
    grid = load_or_build_steric(fw_name, fw_path, max(guest_radius, 1.8))
    n = target_loading(grid, guest_radius)
    res = pack(fw_path, [{"name": guest_name, "file": guest_file, "count": n}],
               use_blocking=True, steric_grid=grid, seed=seed,
               max_iters=max_iters, cutoff_fw=cutoff, cutoff_mm=cutoff,
               blocked_penalty=penalty)
    fw = import_file(fw_path).compute()
    pos_fw = np.array(fw.particles.positions)
    cell = np.array(fw.cell.matrix)[:3, :3]
    inv = np.linalg.inv(cell)
    mols = split_molecules(res.positions, res.n_framework, res.guest_atom_counts)
    m = evaluate(mols, res.guest_centers_frac, pos_fw, cell, inv, grid,
                 cutoff_fw=cutoff, cutoff_mm=cutoff)
    return {
        "framework": fw_name, "guest": guest_name, "loading": n,
        "cutoff_A": cutoff, "blocked_penalty": penalty, "seed": seed,
        "max_iters": max_iters, "converged_iter": res.converged_iter,
        "inaccessible": m["n_inaccessible"], "framework_clash": m["n_gf_clash"],
        "guest_overlap_molecules": m["n_gg_overlap_mols"],
        "total_violations": m["n_inaccessible"] + m["n_gf_clash"] + m["n_gg_overlap_mols"],
        "min_guest_framework_A": float(m["gf_min_dist"].min()),
        "min_guest_guest_A": float(m["gg_nn_dist"].min()) if len(m["gg_nn_dist"]) else None,
    }


def packing_sensitivity():
    affected = (("LTL", "EC"), ("UiO-66", "EC"),
                ("UiO-66", "LiPF6"), ("COF-5", "EC"))
    rows = []
    # Diagnose stochasticity and incomplete convergence at the published cutoff.
    for fw, guest in affected:
        for iters in (12000, 30000):
            for seed in range(5):
                row = run_case(fw, guest, 2.2, 100.0, seed, iters)
                rows.append(row); print("diagnostic", row, flush=True)
    # Cutoff robustness for chemically distinct representative guests and failures.
    cutoff_cases = (("FAU", "CO2"), ("FAU", "H2O"), ("LTL", "EC"),
                    ("UiO-66", "LiPF6"), ("COF-5", "EC"))
    for fw, guest in cutoff_cases:
        for cutoff in (1.8, 2.0, 2.2, 2.4, 2.6):
            for seed in range(3):
                row = run_case(fw, guest, cutoff, 100.0, seed, 30000)
                rows.append(row); print("cutoff", row, flush=True)
    # Penalty sensitivity under a demanding FAU loading.
    for penalty in (0.1, 1.0, 10.0, 100.0, 1000.0):
        for seed in range(5):
            row = run_case("FAU", "CO2", 2.2, penalty, seed, 30000)
            rows.append(row); print("penalty", row, flush=True)
    write_csv("packing_sensitivity.csv", rows)
    return rows


def summarize(grid_rows, pack_rows):
    summary = {"grid_rows": len(grid_rows), "packing_rows": len(pack_rows)}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--section",
        choices=("all", "grid", "packing"),
        default="all",
        help="Run all diagnostics or one independently reproducible section.",
    )
    args = parser.parse_args()
    grids = (
        grid_sensitivity()
        if args.section in ("all", "grid")
        else read_csv("grid_and_filter_sensitivity.csv")
    )
    packs = (
        packing_sensitivity()
        if args.section in ("all", "packing")
        else read_csv("packing_sensitivity.csv")
    )
    summarize(grids, packs)
