#!/usr/bin/env python3
"""Generate probe-radius accessibility sweep data for Figure 2a."""

import numpy as np

from pore_accessibility import load_framework_ovito, probe_sweep

FRAMEWORKS = {
    "FAU": "frameworks/zeolites/FAU.cif",
    "LTL": "frameworks/zeolites/LTL.cif",
    "ERI": "frameworks/zeolites/ERI.cif",
    "OFF": "frameworks/zeolites/OFF.cif",
    "MAZ": "frameworks/zeolites/MAZ.cif",
}

# guest effective radii (kinetic diameter / 2, Angstrom)
GUEST_RADII = {
    "Li+": 0.76, "H2O": 1.33, "CO2": 1.65, "N2": 1.83, "EC": 2.75,
}

RADII = np.round(np.arange(0.6, 3.41, 0.2), 2)
GRID_N = 48
OUT = "runs/probe_sweep.npz"


def main():
    save = {"radii": RADII}
    for name, path in FRAMEWORKS.items():
        print(f"[sweep] {name} ...", flush=True)
        pos, cell, inv_cell, _ = load_framework_ovito(path)
        res = probe_sweep(pos, cell, inv_cell, RADII, grid_n=GRID_N)
        save[f"{name}_acc_void"] = res["accessible_void_fraction"]
        save[f"{name}_acc_cell"] = res["accessible_cell_fraction"]
        save[f"{name}_blocked_void"] = res["blocked_void_fraction"]
        save[f"{name}_largest"] = res["largest_cluster_voxels"]
        save[f"{name}_nclusters"] = res["n_blocked_clusters"]
    save["frameworks"] = np.array(list(FRAMEWORKS.keys()))
    save["guest_names"] = np.array(list(GUEST_RADII.keys()))
    save["guest_radii"] = np.array(list(GUEST_RADII.values()))
    np.savez_compressed(OUT, **save)
    print(f"[done] saved {OUT}")


if __name__ == "__main__":
    main()
