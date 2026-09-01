"""Periodic-cell geometry helpers shared by packing and validation code."""

from __future__ import annotations

from itertools import product

import numpy as np


_NEIGHBOR_SHIFTS = np.asarray(list(product((-1, 0, 1), repeat=3)), dtype=float)


def minimum_image(vecs, cell, inv_cell=None):
    """Return exact shortest Cartesian images for a general 3-D cell.

    The common 27-image search is used first.  A singular-value lower bound
    certifies that result for ordinary reduced cells.  Only uncertified vectors
    fall back to a wider, finite search whose radius is guaranteed to contain
    the closest lattice image.
    """
    vecs = np.asarray(vecs, dtype=float)
    cell = np.asarray(cell, dtype=float)
    if vecs.shape[-1:] != (3,):
        raise ValueError("vecs must have shape (..., 3)")
    if cell.shape != (3, 3):
        raise ValueError("cell must have shape (3, 3)")
    if inv_cell is None:
        inv_cell = np.linalg.inv(cell)
    else:
        inv_cell = np.asarray(inv_cell, dtype=float)

    original_shape = vecs.shape
    flat_frac = (vecs.reshape(-1, 3) @ inv_cell)
    wrapped = flat_frac - np.round(flat_frac)

    best = wrapped @ cell
    best_norm2 = np.einsum("ni,ni->n", best, best)
    for shift in _NEIGHBOR_SHIFTS:
        if not np.any(shift):
            continue
        candidate = (wrapped + shift) @ cell
        candidate_norm2 = np.einsum("ni,ni->n", candidate, candidate)
        improved = candidate_norm2 < best_norm2
        best[improved] = candidate[improved]
        best_norm2[improved] = candidate_norm2[improved]
    best_norm = np.sqrt(best_norm2)

    # Any image outside the 27-image cube has at least one fractional
    # component with magnitude >= 1.5.  The smallest singular value therefore
    # supplies a rigorous Cartesian lower bound for every omitted image.
    sigma_min = float(np.linalg.svd(cell, compute_uv=False)[-1])
    if not np.isfinite(sigma_min) or sigma_min <= 0.0:
        raise ValueError("cell must be finite and nonsingular")
    unresolved = np.flatnonzero(best_norm > (1.5 * sigma_min + 1e-12))

    if unresolved.size:
        # If an image improves on the current best distance U, its fractional
        # norm is at most U/sigma_min.  Thus shifts beyond ceil(U/sigma_min +
        # 0.5) cannot improve the result.  Group by radius to avoid imposing a
        # pathological cell's search size on every vector.
        radii = np.ceil(best_norm[unresolved] / sigma_min + 0.5).astype(int)
        for radius in np.unique(radii):
            indices = unresolved[radii == radius]
            group = wrapped[indices]
            group_best = best[indices].copy()
            group_best_norm2 = np.einsum("ni,ni->n", group_best, group_best)
            for shift in product(range(-radius, radius + 1), repeat=3):
                shift = np.asarray(shift, dtype=float)
                candidate = (group + shift) @ cell
                candidate_norm2 = np.einsum("ni,ni->n", candidate, candidate)
                improved = candidate_norm2 < group_best_norm2
                group_best[improved] = candidate[improved]
                group_best_norm2[improved] = candidate_norm2[improved]
            best[indices] = group_best

    return best.reshape(original_shape)
