#!/usr/bin/env python3
"""Grid-based pore accessibility analysis for periodic frameworks."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_NEIGHBORS_6 = (
    (1, 0, 0), (-1, 0, 0),
    (0, 1, 0), (0, -1, 0),
    (0, 0, 1), (0, 0, -1),
)

_NEIGHBORS_26 = tuple(
    (di, dj, dk)
    for di in (-1, 0, 1)
    for dj in (-1, 0, 1)
    for dk in (-1, 0, 1)
    if (di, dj, dk) != (0, 0, 0)
)


def minimum_image(vecs, cell, inv_cell):
    frac = vecs @ inv_cell
    frac -= np.round(frac)
    return frac @ cell


def min_distance_points_to_framework(cart_pts, pos_fw, cell, inv_cell, batch=2048):
    mins = np.empty(len(cart_pts), dtype=float)
    for start in range(0, len(cart_pts), batch):
        pts = cart_pts[start : start + batch]
        diff = pts[:, None, :] - pos_fw[None, :, :]
        diff = minimum_image(diff, cell, inv_cell)
        mins[start : start + batch] = np.linalg.norm(diff, axis=2).min(axis=1)
    return mins


def _grid_frac_centers(grid_n):
    t = (np.arange(grid_n) + 0.5) / grid_n
    fx, fy, fz = np.meshgrid(t, t, t, indexing="ij")
    return np.stack([fx, fy, fz], axis=-1)


def _apply_manual_blocked(blocked, manual_boxes):
    if not manual_boxes:
        return blocked
    n = blocked.shape[0]
    t = (np.arange(n) + 0.5) / n
    fx, fy, fz = np.meshgrid(t, t, t, indexing="ij")
    for box in manual_boxes:
        fmin = np.asarray(box["frac_min"], dtype=float)
        fmax = np.asarray(box["frac_max"], dtype=float)
        mask = (
            (fx >= fmin[0]) & (fx <= fmax[0])
            & (fy >= fmin[1]) & (fy <= fmax[1])
            & (fz >= fmin[2]) & (fz <= fmax[2])
        )
        blocked |= mask
    return blocked


def _neighbor_steps(connectivity):
    if connectivity == 6:
        return _NEIGHBORS_6
    if connectivity == 26:
        return _NEIGHBORS_26
    raise ValueError("connectivity must be 6 or 26")


def _flood_fill_accessible(void_mask, connectivity=26):
    """BFS flood fill from unit-cell boundary void voxels (periodic)."""
    n = void_mask.shape[0]
    accessible = np.zeros_like(void_mask, dtype=bool)
    steps = _neighbor_steps(connectivity)
    q = []

    def try_add(i, j, k):
        if void_mask[i, j, k] and not accessible[i, j, k]:
            accessible[i, j, k] = True
            q.append((i, j, k))

    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i in (0, n - 1) or j in (0, n - 1) or k in (0, n - 1):
                    try_add(i, j, k)

    while q:
        i, j, k = q.pop()
        for di, dj, dk in steps:
            try_add((i + di) % n, (j + dj) % n, (k + dk) % n)

    return accessible


def _cluster_sizes(mask, connectivity=26):
    n = mask.shape[0]
    steps = _neighbor_steps(connectivity)
    visited = np.zeros_like(mask, dtype=bool)
    sizes = []

    for i in range(n):
        for j in range(n):
            for k in range(n):
                if not mask[i, j, k] or visited[i, j, k]:
                    continue
                stack = [(i, j, k)]
                visited[i, j, k] = True
                count = 0
                while stack:
                    ci, cj, ck = stack.pop()
                    count += 1
                    for di, dj, dk in steps:
                        ni, nj, nk = (ci + di) % n, (cj + dj) % n, (ck + dk) % n
                        if mask[ni, nj, nk] and not visited[ni, nj, nk]:
                            visited[ni, nj, nk] = True
                            stack.append((ni, nj, nk))
                sizes.append(count)
    return sizes


def _filter_small_blocked_clusters(blocked, min_cluster_size, connectivity=26):
    """Drop speckle: tiny blocked clusters are treated as accessible artifacts."""
    if min_cluster_size <= 1:
        return blocked.copy()

    n = blocked.shape[0]
    steps = _neighbor_steps(connectivity)
    visited = np.zeros_like(blocked, dtype=bool)
    filtered = blocked.copy()

    for i in range(n):
        for j in range(n):
            for k in range(n):
                if not blocked[i, j, k] or visited[i, j, k]:
                    continue
                stack = [(i, j, k)]
                component = []
                visited[i, j, k] = True
                while stack:
                    ci, cj, ck = stack.pop()
                    component.append((ci, cj, ck))
                    for di, dj, dk in steps:
                        ni, nj, nk = (ci + di) % n, (cj + dj) % n, (ck + dk) % n
                        if blocked[ni, nj, nk] and not visited[ni, nj, nk]:
                            visited[ni, nj, nk] = True
                            stack.append((ni, nj, nk))
                if len(component) < min_cluster_size:
                    for ci, cj, ck in component:
                        filtered[ci, cj, ck] = False

    return filtered


class PoreGrid:
    """3D grid marking solid / accessible / blocked void voxels."""

    def __init__(
        self,
        blocked_mask,
        accessible_mask,
        solid_mask,
        grid_n,
        probe_radius,
        framework_file="",
        manual_boxes=None,
        connectivity=26,
        min_cluster_size=4,
        mode="steric",
        cluster_stats=None,
    ):
        self.blocked_mask = blocked_mask
        self.accessible_mask = accessible_mask
        self.solid_mask = solid_mask
        self.grid_n = int(grid_n)
        self.probe_radius = float(probe_radius)
        self.framework_file = framework_file
        self.manual_boxes = manual_boxes or []
        self.connectivity = int(connectivity)
        self.min_cluster_size = int(min_cluster_size)
        self.mode = mode
        self.cluster_stats = cluster_stats or {}

        frac = _grid_frac_centers(self.grid_n).reshape(-1, 3)
        acc = accessible_mask.reshape(-1)
        self.accessible_frac_centers = frac[acc]

    @classmethod
    def build(
        cls,
        pos_fw,
        cell,
        inv_cell,
        grid_n=64,
        probe_radius=1.2,
        manual_boxes=None,
        framework_file="",
        connectivity=26,
        min_cluster_size=4,
        mode="topological",
    ):
        frac_grid = _grid_frac_centers(grid_n)
        cart_pts = frac_grid.reshape(-1, 3) @ cell
        dist = min_distance_points_to_framework(cart_pts, pos_fw, cell, inv_cell)
        solid = (dist < probe_radius).reshape(grid_n, grid_n, grid_n)
        void = ~solid
        accessible = _flood_fill_accessible(void, connectivity=connectivity)
        blocked_raw = void & (~accessible)
        blocked = _filter_small_blocked_clusters(
            blocked_raw, min_cluster_size, connectivity=connectivity
        )
        accessible = accessible | (blocked_raw & (~blocked))

        manual_mask = _apply_manual_blocked(np.zeros_like(blocked), manual_boxes)
        blocked = blocked | manual_mask
        accessible = accessible & (~manual_mask)

        sizes = _cluster_sizes(blocked, connectivity=connectivity)
        cluster_stats = {
            "n_clusters": len(sizes),
            "largest_clusters": sorted(sizes, reverse=True)[:10],
            "speckle_removed": int(blocked_raw.sum() - blocked.sum()),
        }

        return cls(
            blocked_mask=blocked,
            accessible_mask=accessible,
            solid_mask=solid,
            grid_n=grid_n,
            probe_radius=probe_radius,
            framework_file=framework_file,
            manual_boxes=manual_boxes or [],
            connectivity=connectivity,
            min_cluster_size=min_cluster_size,
            mode=mode,
            cluster_stats=cluster_stats,
        )

    @property
    def stats(self):
        void = ~self.solid_mask
        return {
            "grid_n": self.grid_n,
            "probe_radius": self.probe_radius,
            "connectivity": self.connectivity,
            "min_cluster_size": self.min_cluster_size,
            "mode": self.mode,
            "solid_voxels": int(self.solid_mask.sum()),
            "accessible_voxels": int(self.accessible_mask.sum()),
            "blocked_voxels": int(self.blocked_mask.sum()),
            "void_voxels": int(void.sum()),
            "accessible_fraction_of_void": float(
                self.accessible_mask.sum() / max(1, void.sum())
            ),
            **self.cluster_stats,
        }

    def print_report(self):
        s = self.stats
        print(f"[PORE] mode={s['mode']}, grid_n={s['grid_n']}, "
              f"probe={s['probe_radius']:.2f} A, conn={s['connectivity']}, "
              f"min_cluster={s['min_cluster_size']}")
        print(f"       solid:      {s['solid_voxels']:8d} / {self.grid_n**3}")
        print(f"       accessible: {s['accessible_voxels']:8d}")
        print(f"       blocked:    {s['blocked_voxels']:8d}")
        print(f"       acc/void:   {100 * s['accessible_fraction_of_void']:.2f}%")
        if s.get("speckle_removed"):
            print(f"       speckle removed: {s['speckle_removed']}")
        if s.get("n_clusters"):
            print(f"       blocked clusters: {s['n_clusters']}, "
                  f"largest={s.get('largest_clusters', [])[:5]}")

    def is_center_blocked(self, frac):
        frac = np.mod(np.asarray(frac, dtype=float), 1.0)
        idx = np.floor(frac * self.grid_n).astype(int) % self.grid_n
        return bool(self.blocked_mask[idx[0], idx[1], idx[2]])

    def is_center_accessible(self, frac):
        frac = np.mod(np.asarray(frac, dtype=float), 1.0)
        idx = np.floor(frac * self.grid_n).astype(int) % self.grid_n
        return bool(self.accessible_mask[idx[0], idx[1], idx[2]])

    def random_accessible_center(self, rng=None):
        rng = np.random if rng is None else rng
        if len(self.accessible_frac_centers) == 0:
            raise RuntimeError("No accessible void voxel available for placement.")
        i = rng.randint(len(self.accessible_frac_centers))
        return self.accessible_frac_centers[i].copy()

    def voxel_centers(self, mask_name="blocked"):
        masks = {
            "blocked": self.blocked_mask,
            "accessible": self.accessible_mask,
            "solid": self.solid_mask,
            "void": ~self.solid_mask,
        }
        frac = _grid_frac_centers(self.grid_n).reshape(-1, 3)
        return frac[masks[mask_name].reshape(-1)]

    def voxel_centers_cart(self, mask_name, cell):
        frac = self.voxel_centers(mask_name)
        return frac @ cell, frac

    def save(self, path):
        path = Path(path)
        meta = {
            "framework_file": self.framework_file,
            "grid_n": self.grid_n,
            "probe_radius": self.probe_radius,
            "connectivity": self.connectivity,
            "min_cluster_size": self.min_cluster_size,
            "mode": self.mode,
            "manual_boxes": self.manual_boxes,
            "stats": self.stats,
        }
        np.savez_compressed(
            path,
            blocked_mask=self.blocked_mask,
            accessible_mask=self.accessible_mask,
            solid_mask=self.solid_mask,
            meta_json=np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
        )

    @classmethod
    def load(cls, path):
        path = Path(path)
        data = np.load(path, allow_pickle=False)
        meta = json.loads(bytes(data["meta_json"].tolist()).decode())
        stats = meta.get("stats", {})
        cluster_stats = {
            k: stats[k]
            for k in ("n_clusters", "largest_clusters", "speckle_removed")
            if k in stats
        }
        return cls(
            blocked_mask=data["blocked_mask"],
            accessible_mask=data["accessible_mask"],
            solid_mask=data["solid_mask"],
            grid_n=meta["grid_n"],
            probe_radius=meta["probe_radius"],
            framework_file=meta.get("framework_file", ""),
            manual_boxes=meta.get("manual_boxes", []),
            connectivity=meta.get("connectivity", 26),
            min_cluster_size=meta.get("min_cluster_size", 4),
            mode=meta.get("mode", "unknown"),
            cluster_stats=cluster_stats,
        )


def load_framework_ovito(path):
    from ovito.io import import_file

    data = import_file(path).compute()
    pos = np.array(data.particles.positions)
    cell = np.array(data.cell.matrix)[:3, :3]
    inv_cell = np.linalg.inv(cell)
    syms = _get_symbols(data)
    return pos, cell, inv_cell, syms


def _get_symbols(data):
    if "Element" in data.particles:
        arr = data.particles["Element"].array
        return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]
    if "Particle Type" in data.particles:
        type_ids = data.particles["Particle Type"].array
        type_map = {t.id: t.name for t in data.particles.particle_types.types}
        return [type_map[i] for i in type_ids]
    raise RuntimeError("No element/type info in framework.")


T_ELEMENTS = frozenset({
    "Si", "Al", "P", "B", "Ga", "Ge", "Ti", "Fe", "Zn", "Be", "T",
})


def build_t_skeleton(pos, syms, cell, inv_cell, t_o_cutoff=2.0):
    """Return T atom coords and T-T edges (two T atoms sharing one O)."""
    t_idx = [i for i, s in enumerate(syms) if s in T_ELEMENTS]
    o_idx = [i for i, s in enumerate(syms) if s == "O"]
    if not t_idx:
        raise RuntimeError("No T atoms (Si/Al/...) found in framework.")

    t_pos = pos[t_idx]
    o_pos = pos[o_idx]
    edges = set()

    for o in o_pos:
        diff = t_pos - o
        diff = minimum_image(diff, cell, inv_cell)
        dist = np.linalg.norm(diff, axis=1)
        order = np.argsort(dist)
        if dist[order[0]] > t_o_cutoff or dist[order[1]] > t_o_cutoff:
            continue
        i, j = int(order[0]), int(order[1])
        if i != j:
            edges.add((min(i, j), max(i, j)))

    return t_pos, sorted(edges)


def t_edge_segment(p1, p2, cell, inv_cell):
    """Cartesian segment between two points with minimum-image convention."""
    frac = (p2 - p1) @ inv_cell
    frac -= np.round(frac)
    return p1, p1 + frac @ cell

