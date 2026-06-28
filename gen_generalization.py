#!/usr/bin/env python3
"""Generate cross-framework generalization data for Figure 4."""

import os

import numpy as np
from ovito.io import import_file

from baselines.run_packmol import packmol_pack
from error_metrics import evaluate, split_molecules
from mc_engine import pack
from pore_accessibility import PoreGrid, load_framework_ovito

FRAMEWORKS = {
    "FAU": ("frameworks/zeolites/FAU.cif", "zeolite"),
    "LTL": ("frameworks/zeolites/LTL.cif", "zeolite"),
    "ERI": ("frameworks/zeolites/ERI.cif", "zeolite"),
    "MOF-5": ("frameworks/mofs/MOF-5.cif", "MOF"),
    "UiO-66": ("frameworks/mofs/UiO-66.cif", "MOF"),
    "COF-5": ("frameworks/cofs/COF-5.cif", "COF"),
}

GUESTS = {
    "CO2": ("molecules/CO2.xyz", 1.65),
    "H2O": ("molecules/H2O.xyz", 1.33),
    "EC": ("molecules/EC.xyz", 2.75),
    "LiPF6": ("molecules/LiPF6.xyz", 3.0),
}

PROBE = 1.8
GRID_N = 48
MAX_ITERS = 12000
CUTOFF = 2.2
OUT = "runs/generalization.npz"

# Packmol only for orthogonal cubic cells
PACKMOL_FW = {"FAU", "MOF-5"}


def is_orthogonal_cubic(cell, tol=1e-2):
    off = np.abs(cell - np.diag(np.diag(cell)))
    return off.max() < tol


def load_or_build_steric(name, framework_file, probe):
    cache = f"runs/steric_{name}_r{probe:.1f}.npz"
    if os.path.isfile(cache):
        return PoreGrid.load(cache)
    pos, cell, inv, _ = load_framework_ovito(framework_file)
    grid = PoreGrid.build(
        pos, cell, inv, grid_n=GRID_N, probe_radius=probe,
        framework_file=framework_file, mode="steric",
    )
    grid.save(cache)
    return grid


def target_loading(grid, guest_radius):
    """Scale loading inversely with guest size relative to CO2."""
    acc = grid.stats["accessible_voxels"]
    scale = (1.65 / guest_radius) ** 2
    n = max(4, int(acc / 700 * scale))
    return min(n, 24)


def main():
    fw_names = list(FRAMEWORKS.keys())
    guest_names = list(GUESTS.keys())
    n_fw, n_g = len(fw_names), len(guest_names)

    loading = np.zeros((n_fw, n_g), dtype=int)
    violations = np.zeros((n_fw, n_g), dtype=int)
    acc_frac = np.zeros(n_fw)
    pk_inacc = np.full((n_fw, n_g), np.nan)
    families = []

    for fi, fw_name in enumerate(fw_names):
        fw_path, family = FRAMEWORKS[fw_name]
        families.append(family)
        grid = load_or_build_steric(fw_name, fw_path, PROBE)
        acc_frac[fi] = grid.stats["accessible_voxels"] / (GRID_N ** 3)

        fw = import_file(fw_path).compute()
        pos_fw = np.array(fw.particles.positions)
        cell = np.array(fw.cell.matrix)[:3, :3]
        inv_cell = np.linalg.inv(cell)
        cubic = is_orthogonal_cubic(cell)

        for gi, gname in enumerate(guest_names):
            gfile, g_radius = GUESTS[gname]
            # guest-size-matched steric grid for fair accessibility judging
            ggrid = load_or_build_steric(fw_name, fw_path, max(g_radius, PROBE))
            n = target_loading(ggrid, g_radius)
            loading[fi, gi] = n

            res = pack(
                fw_path, [{"name": gname, "file": gfile, "count": n}],
                use_blocking=True, steric_grid=ggrid, seed=0,
                max_iters=MAX_ITERS, cutoff_fw=CUTOFF, cutoff_mm=CUTOFF,
            )
            mols = split_molecules(res.positions, res.n_framework, res.guest_atom_counts)
            m = evaluate(mols, res.guest_centers_frac, pos_fw, cell, inv_cell,
                         ggrid, cutoff_fw=CUTOFF, cutoff_mm=CUTOFF)
            violations[fi, gi] = (
                m["n_inaccessible"] + m["n_gf_clash"] + m["n_gg_overlap_mols"]
            )
            print(f"[ours] {fw_name} x {gname}: N={n} viol={violations[fi, gi]}", flush=True)

            if cubic and fw_name in PACKMOL_FW:
                try:
                    pk = packmol_pack(
                        fw_path, [{"name": gname, "file": gfile, "count": n}],
                        cell, tolerance=CUTOFF, seed=42,
                    )
                    pm = evaluate(pk["mol_positions"], pk["centers_frac"],
                                  pos_fw, cell, inv_cell, ggrid,
                                  cutoff_fw=CUTOFF, cutoff_mm=CUTOFF)
                    pk_inacc[fi, gi] = pm["n_inaccessible"]
                    print(f"[packmol] {fw_name} x {gname}: inacc={pm['n_inaccessible']}", flush=True)
                except Exception as exc:
                    print(f"[packmol] {fw_name} x {gname} failed: {exc}", flush=True)

    np.savez_compressed(
        OUT,
        frameworks=np.array(fw_names),
        guests=np.array(guest_names),
        families=np.array(families),
        loading=loading,
        violations=violations,
        accessible_cell_fraction=acc_frac,
        packmol_inaccessible=pk_inacc,
    )
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
