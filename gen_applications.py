#!/usr/bin/env python3
"""Generate Figure 3 application data: our method vs Packmol on 3 cases.

Cases (all in FAU, cubic -> valid for Packmol's box constraint):
  1. Li-salt electrolyte : EC + DMC + LiPF6 at realistic ratio
  2. CO2 adsorption      : CO2 loading
  3. CO2/H2O separation  : competitive CO2 + H2O mixture

For each case we run our accessibility-aware MC and Packmol (accessibility- and
PBC-blind), score both with the same physical-violation judge, and cache the
structures (for OVITO rendering) plus metrics.
"""

import numpy as np
from ovito.io import import_file

from mc_engine import pack
from error_metrics import evaluate, split_molecules
from baselines.run_packmol import packmol_pack
from pore_accessibility import PoreGrid

FRAMEWORK = "frameworks/zeolites/FAU.cif"
STERIC = "frameworks/zeolites/FAU.blocked.npz"

CASES = {
    "electrolyte": [
        {"name": "EC", "file": "molecules/EC.xyz", "count": 8},
        {"name": "DMC", "file": "molecules/DMC.xyz", "count": 8},
        {"name": "LiPF6", "file": "molecules/LiPF6.xyz", "count": 4},
    ],
    "co2": [
        {"name": "CO2", "file": "molecules/CO2.xyz", "count": 32},
    ],
    "separation": [
        {"name": "CO2", "file": "molecules/CO2.xyz", "count": 16},
        {"name": "H2O", "file": "molecules/H2O.xyz", "count": 16},
    ],
}

SEED = 0
MAX_ITERS = 15000
CUTOFF = 2.2


def main():
    fw = import_file(FRAMEWORK).compute()
    cell = np.array(fw.cell.matrix)[:3, :3]
    inv_cell = np.linalg.inv(cell)
    pos_fw = np.array(fw.particles.positions)
    grid = PoreGrid.load(STERIC)

    summary = {}
    for case, species in CASES.items():
        print(f"\n===== case: {case} =====", flush=True)

        # ---- our method ----
        res = pack(FRAMEWORK, species, use_blocking=True, steric_grid=grid,
                   seed=SEED, max_iters=MAX_ITERS, cutoff_fw=CUTOFF, cutoff_mm=CUTOFF)
        ours_mols = split_molecules(res.positions, res.n_framework, res.guest_atom_counts)
        m_ours = evaluate(ours_mols, res.guest_centers_frac, pos_fw, cell, inv_cell,
                          grid, cutoff_fw=CUTOFF, cutoff_mm=CUTOFF)
        print(f"  ours    : inacc={m_ours['n_inaccessible']} "
              f"gf={m_ours['n_gf_clash']} gg={m_ours['n_gg_overlap_mols']} "
              f"conv_iter={res.converged_iter}")

        # ---- packmol ----
        pk = packmol_pack(FRAMEWORK, species, cell, tolerance=CUTOFF, inset=0.0)
        m_pk = evaluate(pk["mol_positions"], pk["centers_frac"], pos_fw, cell, inv_cell,
                        grid, cutoff_fw=CUTOFF, cutoff_mm=CUTOFF)
        print(f"  packmol : inacc={m_pk['n_inaccessible']} "
              f"gf={m_pk['n_gf_clash']} gg={m_pk['n_gg_overlap_mols']}")

        np.savez_compressed(
            f"runs/app_{case}.npz",
            # structures
            ours_positions=res.positions, ours_symbols=np.array(res.symbols),
            ours_nfw=res.n_framework,
            ours_atom_counts=np.array(res.guest_atom_counts),
            ours_bad=m_ours["bad_mask"],
            pk_positions=pk["positions"], pk_symbols=np.array(pk["symbols"]),
            pk_nfw=pk["n_framework"],
            pk_atom_counts=np.array(pk["guest_atom_counts"]),
            pk_bad=m_pk["bad_mask"],
            cell=cell,
            # metrics
            ours_inacc=m_ours["n_inaccessible"], ours_gf=m_ours["n_gf_clash"],
            ours_gg=m_ours["n_gg_overlap_mols"], ours_n=m_ours["n_molecules"],
            pk_inacc=m_pk["n_inaccessible"], pk_gf=m_pk["n_gf_clash"],
            pk_gg=m_pk["n_gg_overlap_mols"], pk_n=m_pk["n_molecules"],
            ours_gf_min=m_ours["gf_min_dist"], ours_gg_nn=m_ours["gg_nn_dist"],
            pk_gf_min=m_pk["gf_min_dist"], pk_gg_nn=m_pk["gg_nn_dist"],
        )
        summary[case] = {
            "ours": (m_ours["n_inaccessible"], m_ours["n_gf_clash"], m_ours["n_gg_overlap_mols"]),
            "packmol": (m_pk["n_inaccessible"], m_pk["n_gf_clash"], m_pk["n_gg_overlap_mols"]),
            "n": m_ours["n_molecules"],
        }
    print("\n[summary]")
    for c, s in summary.items():
        print(f"  {c}: ours={s['ours']} packmol={s['packmol']} N={s['n']}")
    print("[done] runs/app_*.npz")


if __name__ == "__main__":
    main()
