#!/usr/bin/env python3
"""Reusable simulated-annealing Monte Carlo engine for guest packing.

Refactors the logic of insert_collision_multi.py into a callable `pack()` that
returns the final structure plus diagnostics (energy/temperature/acceptance
traces, convergence iteration, and number of guests misplaced into sterically
inaccessible cages). Used by the figure-generation pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from ovito.io import import_file

from periodic_geometry import minimum_image
from pore_accessibility import PoreGrid, count_centers_in_blocked


# ---------- geometry helpers ----------
def random_unit_vector(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def random_small_rotation(max_angle_deg, rng):
    theta = np.deg2rad(max_angle_deg) * (2 * rng.random() - 1.0)
    axis = random_unit_vector(rng)
    kx, ky, kz = axis
    k = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    return np.eye(3) + np.sin(theta) * k + (1 - np.cos(theta)) * (k @ k)


def random_rotation(rng):
    u1, u2, u3 = rng.random(3)
    q = np.array([
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ])
    w, x, y, z = q
    return np.array([
        [1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * z * w, 2 * x * z + 2 * y * w],
        [2 * x * y + 2 * z * w, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * x * w],
        [2 * x * z - 2 * y * w, 2 * y * z + 2 * x * w, 1 - 2 * x * x - 2 * y * y],
    ])


def get_symbols(data):
    if "Element" in data.particles:
        arr = data.particles["Element"].array
        return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]
    if "Particle Type" in data.particles:
        type_ids = data.particles["Particle Type"].array
        type_map = {t.id: t.name for t in data.particles.particle_types.types}
        return [type_map[i] for i in type_ids]
    raise RuntimeError("No element/type info found in data.")


def min_distance(a, b, cell, inv_cell):
    diff = a[:, None, :] - b[None, :, :]
    diff = minimum_image(diff, cell, inv_cell)
    return float(np.min(np.linalg.norm(diff, axis=2)))


def wrap_to_cell(pos, cell, inv_cell):
    frac = pos @ inv_cell
    frac -= np.floor(frac)
    return frac @ cell


def distribute_by_ratio(ratios, n_total):
    ratios = np.asarray(ratios, dtype=float)
    if np.any(ratios <= 0):
        raise ValueError("All ratios must be positive.")
    ideal = ratios / ratios.sum() * n_total
    counts = np.floor(ideal).astype(int)
    remainder = n_total - counts.sum()
    if remainder:
        order = np.argsort(-(ideal - counts))
        for k in range(remainder):
            counts[order[k % len(counts)]] += 1
    return counts


# ---------- result container ----------
@dataclass
class PackResult:
    positions: np.ndarray          # framework + guests, cartesian
    symbols: list
    cell: np.ndarray
    n_framework: int
    guest_type_ids: np.ndarray     # per-guest species index
    guest_centers_frac: np.ndarray
    species_names: list
    energy_trace: list = field(default_factory=list)
    temperature_trace: list = field(default_factory=list)
    acceptance_trace: list = field(default_factory=list)
    converged_iter: int = -1
    n_misplaced: int = 0
    misplaced_mask: np.ndarray = None
    guest_atom_counts: list = field(default_factory=list)


def _load_species(species_cfg, n_total, default_fw, default_mm):
    has_count = all("count" in sp for sp in species_cfg)
    has_ratio = all("ratio" in sp for sp in species_cfg)
    if has_count:
        counts = [int(sp["count"]) for sp in species_cfg]
    elif has_ratio:
        counts = distribute_by_ratio([sp["ratio"] for sp in species_cfg], n_total).tolist()
    else:
        raise ValueError("Each species must define 'ratio' or 'count'.")

    loaded, type_ids = [], []
    for sp, count in zip(species_cfg, counts):
        if count <= 0:
            continue
        mol = import_file(sp["file"]).compute()
        pos = np.array(mol.particles.positions)
        pos -= pos.mean(axis=0)
        loaded.append({
            "name": sp.get("name", sp["file"]),
            "pos": pos,
            "sym": get_symbols(mol),
            "cutoff_fw": sp.get("cutoff_fw", default_fw),
            "cutoff_mm": sp.get("cutoff_mm", default_mm),
            "count": count,
        })
        type_ids.extend([len(loaded) - 1] * count)
    return loaded, np.asarray(type_ids, dtype=int)


def pack(
    framework_file,
    species,
    n_total=None,
    use_blocking=True,
    blocked_file=None,
    steric_grid=None,
    seed=0,
    max_iters=30000,
    cutoff_fw=2.2,
    cutoff_mm=2.2,
    blocked_penalty=100.0,
    t0=1.0,
    t_min=0.02,
    step_transl=0.8,
    step_rot_deg=20.0,
    p_small=0.70,
    p_jump=0.10,
    e_tol=1e-6,
    record_every=50,
    verbose=False,
):
    """Pack guest molecules into a framework via simulated-annealing MC.

    steric_grid : optional PoreGrid used ONLY to measure misplacement (always
                  evaluated, regardless of use_blocking). If None, falls back to
                  blocked_file or no measurement.
    """
    rng = np.random.default_rng(seed)

    fw = import_file(framework_file).compute()
    pos_fw = np.array(fw.particles.positions)
    sym_fw = get_symbols(fw)
    cell = np.array(fw.cell.matrix)[:3, :3]
    inv_cell = np.linalg.inv(cell)

    # grids: planning grid (for placement/penalty) + measurement grid (always)
    plan_grid = None
    if use_blocking:
        if steric_grid is not None:
            plan_grid = steric_grid
        elif blocked_file:
            plan_grid = PoreGrid.load(blocked_file)
    measure_grid = steric_grid
    if measure_grid is None and blocked_file:
        measure_grid = PoreGrid.load(blocked_file)

    sp_list, type_ids = _load_species(species, n_total, cutoff_fw, cutoff_mm)
    n_mol = len(type_ids)
    cutoffs_fw = np.array([sp_list[t]["cutoff_fw"] for t in type_ids])
    cutoffs_mm = np.array([sp_list[t]["cutoff_mm"] for t in type_ids])

    # framework KD-tree over 3x3x3 images for fast minimum-image guest-framework distance
    try:
        from scipy.spatial import cKDTree
        offsets = np.array([
            i * cell[0] + j * cell[1] + k * cell[2]
            for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)
        ])
        fw_rep = (pos_fw[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
        fw_tree = cKDTree(fw_rep)
    except Exception:
        fw_tree = None

    def rand_center():
        if plan_grid is None:
            return rng.random(3)
        centers = plan_grid.accessible_frac_centers
        return centers[rng.integers(len(centers))].copy()

    centers_frac = np.array([rand_center() for _ in range(n_mol)])
    rot_mats = np.array([random_rotation(rng) for _ in range(n_mol)])

    def mol_world(i):
        return sp_list[type_ids[i]]["pos"] @ rot_mats[i].T + (centers_frac[i] @ cell)

    def fw_dist(world_i):
        if fw_tree is not None:
            return float(fw_tree.query(world_i, k=1)[0].min())
        return min_distance(world_i, pos_fw, cell, inv_cell)

    def e_fw_of(world_i, i):
        d = fw_dist(world_i)
        c = cutoffs_fw[i]
        return (c - d) ** 2 if d < c else 0.0

    def e_block_of(i):
        if plan_grid is not None and plan_grid.is_center_blocked(centers_frac[i]):
            return blocked_penalty
        return 0.0

    def e_pair_of(world_i, world_j, i, j):
        d = min_distance(world_i, world_j, cell, inv_cell)
        c = max(cutoffs_mm[i], cutoffs_mm[j])
        return (c - d) ** 2 if d < c else 0.0

    # initial bookkeeping
    world = [mol_world(i) for i in range(n_mol)]
    e_fw = np.array([e_fw_of(world[i], i) for i in range(n_mol)])
    e_block = np.array([e_block_of(i) for i in range(n_mol)])
    Epair = np.zeros((n_mol, n_mol))
    for i in range(n_mol):
        for j in range(i + 1, n_mol):
            v = e_pair_of(world[i], world[j], i, j)
            Epair[i, j] = Epair[j, i] = v

    def conflicts_vec():
        return e_fw + e_block + Epair.sum(axis=1)

    conflicts = conflicts_vec()
    e = float(e_fw.sum() + e_block.sum() + np.triu(Epair, 1).sum())

    e_trace, t_trace, acc_trace = [], [], []
    n_acc = 0
    converged = -1

    for it in range(1, max_iters + 1):
        t = max(t_min, t0 * (1.0 - it / max_iters))
        if e < e_tol:
            converged = it
            break

        w = int(np.argmax(conflicts))
        frac_old = centers_frac[w].copy()
        r_old = rot_mats[w].copy()

        r = rng.random()
        if r < p_small:
            d_cart = step_transl * (2 * rng.random() - 1) * random_unit_vector(rng)
            centers_frac[w] = (frac_old + d_cart @ inv_cell) % 1.0
            rot_mats[w] = random_small_rotation(step_rot_deg, rng) @ r_old
        elif r < p_small + p_jump:
            centers_frac[w] = rand_center()
            rot_mats[w] = random_rotation(rng)
        else:
            centers_frac[w] = frac_old
            rot_mats[w] = random_small_rotation(60.0, rng) @ r_old

        # incremental energy of moved molecule only
        world_w = mol_world(w)
        new_e_fw_w = e_fw_of(world_w, w)
        new_e_block_w = e_block_of(w)
        new_pair_w = np.zeros(n_mol)
        for j in range(n_mol):
            if j == w:
                continue
            new_pair_w[j] = e_pair_of(world_w, world[j], w, j)

        de = ((new_e_fw_w - e_fw[w]) + (new_e_block_w - e_block[w])
              + (new_pair_w.sum() - Epair[w].sum()))

        accept = de <= 0.0 or rng.random() < np.exp(-de / t)
        if accept:
            world[w] = world_w
            e_fw[w] = new_e_fw_w
            e_block[w] = new_e_block_w
            Epair[w, :] = new_pair_w
            Epair[:, w] = new_pair_w
            e += de
            conflicts = conflicts_vec()
            n_acc += 1
        else:
            centers_frac[w] = frac_old
            rot_mats[w] = r_old

        if it % record_every == 0:
            e_trace.append(max(e, 0.0))
            t_trace.append(t)
            acc_trace.append(n_acc / it)
            if verbose:
                print(f"[{it}] E={e:.4f} T={t:.3f} acc={n_acc/it:.2f}")

    # assemble structure
    all_pos = pos_fw.copy()
    all_sym = list(sym_fw)
    for i in range(n_mol):
        all_pos = np.vstack((all_pos, world[i]))
        all_sym.extend(sp_list[type_ids[i]]["sym"])
    all_pos = wrap_to_cell(all_pos, cell, inv_cell)

    # misplacement measured against steric grid (always, even if blocking off)
    n_mis, mis_mask = 0, np.zeros(n_mol, dtype=bool)
    if measure_grid is not None:
        n_mis, mis_mask = count_centers_in_blocked(centers_frac, measure_grid)

    return PackResult(
        positions=all_pos,
        symbols=all_sym,
        cell=cell,
        n_framework=len(pos_fw),
        guest_type_ids=type_ids,
        guest_centers_frac=centers_frac,
        species_names=[s["name"] for s in sp_list],
        energy_trace=e_trace,
        temperature_trace=t_trace,
        acceptance_trace=acc_trace,
        converged_iter=converged,
        n_misplaced=n_mis,
        misplaced_mask=mis_mask,
        guest_atom_counts=[len(sp_list[type_ids[i]]["pos"]) for i in range(n_mol)],
    )
