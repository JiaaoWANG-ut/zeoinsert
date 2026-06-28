#!/usr/bin/env python3
"""Unified physical-violation metrics, shared by all packing methods.

Given a set of guest molecules placed in a periodic framework, quantify:
  1. inaccessible : guest centers landing in sterically blocked cages
  2. gg_overlap   : guest-guest pairs closer than cutoff_mm (PBC)
  3. gf_clash     : guests penetrating the framework (< cutoff_fw, PBC)

The same judge is applied to our method, Packmol, and random baselines so the
comparison in Figure 3 is fair.
"""

from __future__ import annotations

import numpy as np

from pore_accessibility import count_centers_in_blocked


def minimum_image(vecs, cell, inv_cell):
    frac = vecs @ inv_cell
    frac -= np.round(frac)
    return frac @ cell


def _min_dist(a, b, cell, inv_cell):
    diff = a[:, None, :] - b[None, :, :]
    diff = minimum_image(diff, cell, inv_cell)
    return float(np.min(np.linalg.norm(diff, axis=2)))


def evaluate(
    mol_positions,
    centers_frac,
    pos_fw,
    cell,
    inv_cell,
    grid,
    cutoff_fw=2.2,
    cutoff_mm=2.2,
):
    """Evaluate violations.

    mol_positions : list of (n_atoms_i, 3) cartesian arrays, one per molecule.
    centers_frac  : (N, 3) fractional centers of molecules.
    grid          : steric PoreGrid for accessibility judging (or None).
    Returns a dict of counts, fractions, masks, and min-distance arrays.
    """
    n = len(mol_positions)

    # 1. inaccessible cages
    if grid is not None and n:
        n_inacc, inacc_mask = count_centers_in_blocked(centers_frac, grid)
    else:
        n_inacc, inacc_mask = 0, np.zeros(n, dtype=bool)

    # 3. guest-framework clash + min distances
    gf_min = np.array([_min_dist(m, pos_fw, cell, inv_cell) for m in mol_positions]) \
        if n else np.array([])
    gf_mask = gf_min < cutoff_fw
    n_gf = int(gf_mask.sum())

    # 2. guest-guest overlap + nearest-neighbour distances
    gg_pair_count = 0
    gg_mol_mask = np.zeros(n, dtype=bool)
    gg_nn = np.full(n, np.inf)
    for i in range(n):
        for j in range(i + 1, n):
            d = _min_dist(mol_positions[i], mol_positions[j], cell, inv_cell)
            gg_nn[i] = min(gg_nn[i], d)
            gg_nn[j] = min(gg_nn[j], d)
            if d < cutoff_mm:
                gg_pair_count += 1
                gg_mol_mask[i] = True
                gg_mol_mask[j] = True
    gg_nn = gg_nn[np.isfinite(gg_nn)]

    # any molecule that violates anything (for render highlighting)
    bad_mask = inacc_mask | gf_mask | gg_mol_mask

    denom = max(1, n)
    return {
        "n_molecules": n,
        "n_inaccessible": n_inacc,
        "frac_inaccessible": n_inacc / denom,
        "n_gf_clash": n_gf,
        "frac_gf_clash": n_gf / denom,
        "n_gg_overlap_pairs": gg_pair_count,
        "n_gg_overlap_mols": int(gg_mol_mask.sum()),
        "frac_gg_overlap": int(gg_mol_mask.sum()) / denom,
        "inaccessible_mask": inacc_mask,
        "gf_clash_mask": gf_mask,
        "gg_overlap_mask": gg_mol_mask,
        "bad_mask": bad_mask,
        "gf_min_dist": gf_min,
        "gg_nn_dist": gg_nn,
    }


def split_molecules(positions, n_framework, atom_counts):
    """Split a combined (framework+guests) cartesian array into per-molecule lists."""
    guests = positions[n_framework:]
    mols, start = [], 0
    for c in atom_counts:
        mols.append(guests[start:start + c])
        start += c
    return mols
