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
    """Nearest framework-atom distance for each point (minimum-image)."""
    try:
        from scipy.spatial import cKDTree
        # replicate framework into 3x3x3 images so a plain KD-tree captures
        # the minimum-image distance for any point inside the cell
        offsets = np.array([
            i * cell[0] + j * cell[1] + k * cell[2]
            for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)
        ])
        rep = (pos_fw[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
        tree = cKDTree(rep)
        dist, _ = tree.query(cart_pts, k=1, workers=-1)
        return dist
    except Exception:
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


def _structure_element(connectivity):
    from scipy import ndimage
    return ndimage.generate_binary_structure(3, 1 if connectivity == 6 else 3)


def _connected_components(mask, connectivity=26):
    """Label periodic 3D connected components using scipy + boundary merge.

    scipy.ndimage.label is non-periodic; we merge labels that touch across
    opposite faces via union-find to recover periodicity.
    """
    from scipy import ndimage

    struct = _structure_element(connectivity)
    labels, n = ndimage.label(mask, structure=struct)
    if n == 0:
        return labels.astype(np.int32), []

    parent = list(range(n + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    # merge labels adjacent across periodic neighbors (unique pairs only)
    shifts = _NEIGHBORS_6 if connectivity == 6 else _NEIGHBORS_26
    for (di, dj, dk) in shifts:
        rolled = np.roll(labels, (-di, -dj, -dk), axis=(0, 1, 2))
        touch = (labels > 0) & (rolled > 0) & (labels != rolled)
        if not touch.any():
            continue
        pairs = np.unique(np.stack([labels[touch], rolled[touch]], axis=1), axis=0)
        for x, y in pairs.tolist():
            union(int(x), int(y))

    # relabel by root
    remap = {}
    new_labels = np.zeros_like(labels, dtype=np.int32)
    nz = labels > 0
    roots = np.vectorize(find)(labels[nz])
    next_id = 0
    sizes_map = {}
    out_ids = np.empty(roots.shape, dtype=np.int32)
    for idx, r in enumerate(roots.tolist()):
        if r not in remap:
            next_id += 1
            remap[r] = next_id
        cid = remap[r]
        out_ids[idx] = cid
        sizes_map[cid] = sizes_map.get(cid, 0) + 1
    new_labels[nz] = out_ids
    sizes = [sizes_map[i] for i in range(1, next_id + 1)]
    return new_labels, sizes


def _component_touches_face(labels, comp_id):
    faces = (
        labels[0, :, :], labels[-1, :, :],
        labels[:, 0, :], labels[:, -1, :],
        labels[:, :, 0], labels[:, :, -1],
    )
    return any((f == comp_id).any() for f in faces)


def _flood_fill_accessible(void_mask, connectivity=26):
    """Accessible void = periodic void components that reach the cell boundary."""
    labels, sizes = _connected_components(void_mask, connectivity=connectivity)
    accessible = np.zeros_like(void_mask, dtype=bool)
    for cid in range(1, len(sizes) + 1):
        if _component_touches_face(labels, cid):
            accessible |= labels == cid
    return accessible


def _cluster_sizes(mask, connectivity=26):
    return _connected_components(mask, connectivity=connectivity)[1]


def _filter_small_blocked_clusters(blocked, min_cluster_size, connectivity=26):
    """Drop speckle: tiny blocked clusters are treated as accessible artifacts."""
    if min_cluster_size <= 1:
        return blocked.copy()
    labels, sizes = _connected_components(blocked, connectivity=connectivity)
    filtered = blocked.copy()
    for lab, sz in enumerate(sizes, start=1):
        if sz < min_cluster_size:
            filtered[labels == lab] = False
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

    def blocked_blobs(self, cell):
        """Collapse each blocked cluster into one (center_cart, radius) blob."""
        labels, sizes = _connected_components(self.blocked_mask,
                                              connectivity=self.connectivity)
        if not sizes:
            return np.zeros((0, 3)), np.zeros(0)
        t = (np.arange(self.grid_n) + 0.5) / self.grid_n
        fx, fy, fz = np.meshgrid(t, t, t, indexing="ij")
        voxel = max(np.linalg.norm(cell, axis=1)) / self.grid_n
        centers, radii = [], []
        for cid, sz in enumerate(sizes, start=1):
            m = labels == cid
            # average via complex exponent to respect periodicity per axis
            cf = []
            for fa in (fx, fy, fz):
                ang = 2 * np.pi * fa[m]
                cf.append((np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())
                           / (2 * np.pi)) % 1.0)
            centers.append(np.array(cf) @ cell)
            radii.append((3 * sz / (4 * np.pi)) ** (1 / 3) * voxel)
        return np.array(centers), np.array(radii)

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


def probe_sweep(pos_fw, cell, inv_cell, radii, grid_n=48,
                connectivity=26, min_cluster_size=4):
    """Scan accessibility as a function of probe radius.

    Returns dict of arrays aligned to `radii`:
      accessible_void_fraction : accessible / void voxels
      accessible_cell_fraction : accessible / total voxels (usable volume)
      blocked_void_fraction    : blocked / void voxels
      largest_cluster_voxels   : size of the largest blocked cluster
      n_blocked_clusters       : number of blocked clusters
    """
    radii = np.asarray(radii, dtype=float)
    total = grid_n ** 3
    out = {
        "radii": radii,
        "accessible_void_fraction": np.zeros_like(radii),
        "accessible_cell_fraction": np.zeros_like(radii),
        "blocked_void_fraction": np.zeros_like(radii),
        "largest_cluster_voxels": np.zeros_like(radii),
        "n_blocked_clusters": np.zeros_like(radii),
    }
    for k, r in enumerate(radii):
        grid = PoreGrid.build(
            pos_fw, cell, inv_cell,
            grid_n=grid_n, probe_radius=float(r),
            connectivity=connectivity, min_cluster_size=min_cluster_size,
            mode="sweep",
        )
        s = grid.stats
        void = max(1, s["void_voxels"])
        largest = s.get("largest_clusters", [])
        out["accessible_void_fraction"][k] = s["accessible_voxels"] / void
        out["accessible_cell_fraction"][k] = s["accessible_voxels"] / total
        out["blocked_void_fraction"][k] = s["blocked_voxels"] / void
        out["largest_cluster_voxels"][k] = largest[0] if largest else 0
        out["n_blocked_clusters"][k] = s.get("n_clusters", 0)
    return out


def count_centers_in_blocked(centers_frac, grid):
    """Count molecule centers (fractional coords, (N,3)) landing in blocked voxels."""
    centers_frac = np.atleast_2d(np.asarray(centers_frac, dtype=float))
    n = grid.grid_n
    idx = np.floor(np.mod(centers_frac, 1.0) * n).astype(int) % n
    flags = grid.blocked_mask[idx[:, 0], idx[:, 1], idx[:, 2]]
    return int(np.count_nonzero(flags)), flags.astype(bool)

