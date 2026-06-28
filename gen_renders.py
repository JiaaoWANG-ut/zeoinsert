#!/usr/bin/env python3
"""Render all OVITO structure panels (PNG) used by Figures 1-4."""

import os
import sys

import numpy as np
from ovito.io import import_file

from mc_engine import pack, get_symbols
from render_ovito import render_structure
from pore_accessibility import PoreGrid, load_framework_ovito

FRAMEWORK = "frameworks/zeolites/FAU.cif"
STERIC = "frameworks/zeolites/FAU.blocked.npz"
PANELS = "figures/panels"
CAM = (-1.0, -1.15, -0.78)
SIZE = (1000, 1000)
FORCE = "--force" in sys.argv

HERO_CASES = [
    ("MOF-5", "frameworks/mofs/MOF-5.cif", "runs/steric_MOF-5_r1.8.npz", "gen_mof5.png"),
    ("UiO-66", "frameworks/mofs/UiO-66.cif", "runs/steric_UiO-66_r1.8.npz", "gen_uio66.png"),
    ("COF-5", "frameworks/cofs/COF-5.cif", "runs/steric_COF-5_r1.8.npz", "gen_cof5.png"),
]


def _skip(path):
    return not FORCE and os.path.isfile(path) and os.path.getsize(path) > 1000


def load_fw(path=FRAMEWORK):
    fw = import_file(path).compute()
    pos = np.array(fw.particles.positions)
    sym = get_symbols(fw)
    cell = np.array(fw.cell.matrix)[:3, :3]
    return pos, sym, cell


def _load_or_build_steric(path, framework_file, probe=1.8):
    if os.path.isfile(path) and not FORCE:
        return PoreGrid.load(path)
    pos, cell, inv, _ = load_framework_ovito(framework_file)
    grid = PoreGrid.build(pos, cell, inv, grid_n=48, probe_radius=probe,
                          framework_file=framework_file, mode="steric")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    grid.save(path)
    return grid


def main():
    pos_fw, sym_fw, cell = load_fw()
    nfw = len(pos_fw)
    grid = PoreGrid.load(STERIC)
    blob_centers, blob_radii = grid.blocked_blobs(cell)
    print(f"[blobs] {len(blob_centers)} blocked cages, radii~{np.round(blob_radii,2)}")

    def render(path, **kwargs):
        if _skip(path):
            print(f"[skip] {path}")
            return
        print(f"[render] {path} ...", flush=True)
        render_structure(**kwargs, out_png=path, size=SIZE, camera_dir=CAM)

    # ---- Fig1: framework only (workflow start) ----
    render(f"{PANELS}/fw_only.png",
           positions=pos_fw, symbols=sym_fw, cell=cell, n_framework=nfw)

    # ---- Fig1: steric accessibility (blocked beta-cages as red blobs) ----
    render(f"{PANELS}/fw_steric_blocked.png",
           positions=pos_fw, symbols=sym_fw, cell=cell, n_framework=nfw,
           blocked_cart=blob_centers, blocked_radii=blob_radii,
           blocked_transparency=0.45)

    # ---- Fig1 + Fig2c: packed result (blocking ON) ----
    res_on = pack(FRAMEWORK, [{"name": "CO2", "file": "molecules/CO2.xyz", "count": 32}],
                  use_blocking=True, steric_grid=grid, seed=0, max_iters=15000)
    render(f"{PANELS}/packed_on.png",
           positions=res_on.positions, symbols=res_on.symbols, cell=res_on.cell,
           n_framework=res_on.n_framework,
           guest_atom_counts=res_on.guest_atom_counts)

    # ---- Fig1 step 4: validated model (guests + blocked cages overlay) ----
    render(f"{PANELS}/packed_validated.png",
           positions=res_on.positions, symbols=res_on.symbols, cell=res_on.cell,
           n_framework=res_on.n_framework,
           guest_atom_counts=res_on.guest_atom_counts,
           blocked_cart=blob_centers, blocked_radii=blob_radii,
           blocked_transparency=0.35)

    # ---- Fig2c: blocking OFF (misplaced guests highlighted) ----
    best = None
    for seed in range(8):
        r = pack(FRAMEWORK, [{"name": "CO2", "file": "molecules/CO2.xyz", "count": 32}],
                 use_blocking=False, steric_grid=grid, seed=seed, max_iters=8000)
        if best is None or r.n_misplaced > best.n_misplaced:
            best = r
        if r.n_misplaced >= 3:
            break
    print(f"[fig2c] blocking OFF misplaced={best.n_misplaced}", flush=True)
    render(f"{PANELS}/packed_off.png",
           positions=best.positions, symbols=best.symbols, cell=best.cell,
           n_framework=best.n_framework,
           guest_atom_counts=best.guest_atom_counts,
           bad_mask=best.misplaced_mask)

    # ---- Fig3: applications, ours vs packmol ----
    for case in ("electrolyte", "co2", "separation"):
        d = np.load(f"runs/app_{case}.npz", allow_pickle=True)
        render(f"{PANELS}/app_{case}_ours.png",
               positions=d["ours_positions"], symbols=list(d["ours_symbols"]),
               cell=d["cell"], n_framework=int(d["ours_nfw"]),
               guest_atom_counts=list(d["ours_atom_counts"]),
               bad_mask=d["ours_bad"])
        render(f"{PANELS}/app_{case}_packmol.png",
               positions=d["pk_positions"], symbols=list(d["pk_symbols"]),
               cell=d["cell"], n_framework=int(d["pk_nfw"]),
               guest_atom_counts=list(d["pk_atom_counts"]),
               bad_mask=d["pk_bad"])
        print(f"[fig3] done {case}", flush=True)

    # ---- Fig4: cross-framework hero renders ----
    for name, fw_path, steric_path, png_name in HERO_CASES:
        if not os.path.isfile(steric_path):
            print(f"[warn] missing {steric_path}, run gen_generalization.py first")
            continue
        sg = PoreGrid.load(steric_path)
        n_target = max(4, int(sg.stats["accessible_voxels"] / 800))
        n_target = min(n_target, 24)
        res = pack(fw_path, [{"name": "CO2", "file": "molecules/CO2.xyz", "count": n_target}],
                   use_blocking=True, steric_grid=sg, seed=0, max_iters=12000)
        bond = min(2.2, max(res.cell.diagonal()) * 0.04) if name == "COF-5" else 2.0
        render(f"{PANELS}/{png_name}",
               positions=res.positions, symbols=res.symbols, cell=res.cell,
               n_framework=res.n_framework,
               guest_atom_counts=res.guest_atom_counts,
               bond_cutoff=bond)

    print("[done] panels in", PANELS)


if __name__ == "__main__":
    main()
