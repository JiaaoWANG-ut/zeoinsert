#!/usr/bin/env python3
"""Visualize blocked / accessible pore regions (topological + steric maps)."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from pore_accessibility import (
    PoreGrid,
    build_t_skeleton,
    load_framework_ovito,
    t_edge_segment,
)

FRAMEWORK_FILE = "frameworks/zeolites/FAU.cif"
BLOCKED_TOPO = "frameworks/zeolites/FAU.blocked_topo.npz"
BLOCKED_STERIC = "frameworks/zeolites/FAU.blocked.npz"

OUTPUT_PNG = "viz/FAU_pores.png"
OUTPUT_CIF_TOPO = "viz/FAU_pores_topo.cif"
OUTPUT_CIF_STERIC = "viz/FAU_pores_steric.cif"

SHOW_FRAMEWORK = False
SHOW_T_SKELETON = True
T_NODE_SIZE = 8
T_EDGE_COLOR = "#2ca02c"
T_EDGE_ALPHA = 0.55
T_EDGE_LW = 0.6
ACCESSIBLE_SAMPLE = 400


def cell_angles(cell):
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = np.degrees(np.arccos(np.dot(cell[1], cell[2]) / (b * c)))
    beta = np.degrees(np.arccos(np.dot(cell[0], cell[2]) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(cell[0], cell[1]) / (a * b)))
    return a, b, c, alpha, beta, gamma


def write_cif(path, cell, inv_cell, symbols, positions):
    a, b, c, alpha, beta, gamma = cell_angles(cell)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("data_pore_viz\n_symmetry_space_group_name_H-M 'P1'\n_symmetry_Int_Tables_number 1\n\n")
        f.write(f"_cell_length_a {a:.6f}\n_cell_length_b {b:.6f}\n_cell_length_c {c:.6f}\n")
        f.write(f"_cell_angle_alpha {alpha:.6f}\n_cell_angle_beta {beta:.6f}\n")
        f.write(f"_cell_angle_gamma {gamma:.6f}\n\nloop_\n")
        f.write("_atom_site_label\n_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        for idx, (sym, pos) in enumerate(zip(symbols, positions), start=1):
            frac = pos @ inv_cell
            f.write(f"{sym}{idx} {sym} {frac[0]:.6f} {frac[1]:.6f} {frac[2]:.6f}\n")


def draw_cell_box(ax, cell):
    o = np.zeros(3)
    for i in range(3):
        e = cell[i]
        ax.plot([o[0], e[0]], [o[1], e[1]], [o[2], e[2]], "k-", lw=0.6, alpha=0.4)


def draw_t_skeleton(ax, t_pos, edges, cell, inv_cell):
    xs, ys, zs = [], [], []
    for i, j in edges:
        p1, p2 = t_edge_segment(t_pos[i], t_pos[j], cell, inv_cell)
        xs.extend([p1[0], p2[0], np.nan])
        ys.extend([p1[1], p2[1], np.nan])
        zs.extend([p1[2], p2[2], np.nan])
    ax.plot(
        xs, ys, zs,
        color=T_EDGE_COLOR, lw=T_EDGE_LW, alpha=T_EDGE_ALPHA,
        label=f"T skeleton ({len(edges)} edges)",
    )
    if T_NODE_SIZE > 0:
        ax.scatter(
            t_pos[:, 0], t_pos[:, 1], t_pos[:, 2],
            c=T_EDGE_COLOR, s=T_NODE_SIZE, alpha=0.9, depthshade=False,
        )


def plot_panel(ax, cell, inv_cell, blocked_cart, title, t_pos=None, t_edges=None,
               accessible_cart=None):
    if SHOW_T_SKELETON and t_pos is not None and t_edges is not None:
        draw_t_skeleton(ax, t_pos, t_edges, cell, inv_cell)

    if len(blocked_cart):
        ax.scatter(
            blocked_cart[:, 0], blocked_cart[:, 1], blocked_cart[:, 2],
            c="red", s=22, alpha=0.9, depthshade=True,
            label=f"blocked ({len(blocked_cart)})",
        )
    if accessible_cart is not None and len(accessible_cart):
        ax.scatter(
            accessible_cart[:, 0], accessible_cart[:, 1], accessible_cart[:, 2],
            c="royalblue", s=3, alpha=0.15,
            label=f"accessible sample ({len(accessible_cart)})",
        )
    draw_cell_box(ax, cell)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)


def export_markers(grid, cell, inv_cell, blocked_sym, out_cif, fw_path):
    blocked_cart, _ = grid.voxel_centers_cart("blocked", cell)
    symbols = []
    positions = []
    if SHOW_FRAMEWORK:
        pos, _, _, syms = load_framework_ovito(fw_path)
        symbols.extend(syms)
        positions.append(pos)
    if len(blocked_cart):
        symbols += [blocked_sym] * len(blocked_cart)
        positions.append(blocked_cart)
    if positions:
        write_cif(out_cif, cell, inv_cell, symbols, np.vstack(positions))
    return blocked_cart


def main():
    pos_fw, cell, inv_cell, syms = load_framework_ovito(FRAMEWORK_FILE)
    t_pos, t_edges = build_t_skeleton(pos_fw, syms, cell, inv_cell)
    print(f"[INFO] T skeleton: {len(t_pos)} T atoms, {len(t_edges)} T-T edges")

    print("[INFO] Loading pore maps...")
    grid_topo = PoreGrid.load(BLOCKED_TOPO)
    grid_steric = PoreGrid.load(BLOCKED_STERIC)
    grid_topo.print_report()
    grid_steric.print_report()

    blocked_topo = export_markers(grid_topo, cell, inv_cell, "Bk", OUTPUT_CIF_TOPO, FRAMEWORK_FILE)
    blocked_steric = export_markers(grid_steric, cell, inv_cell, "Bs", OUTPUT_CIF_STERIC, FRAMEWORK_FILE)
    print(f"[INFO] Saved {OUTPUT_CIF_TOPO} ({len(blocked_topo)} Bk markers)")
    print(f"[INFO] Saved {OUTPUT_CIF_STERIC} ({len(blocked_steric)} Bs markers)")

    acc_frac = grid_steric.voxel_centers("accessible")
    rng = np.random.default_rng(0)
    n = min(ACCESSIBLE_SAMPLE, len(acc_frac))
    acc_sample = acc_frac[rng.choice(len(acc_frac), size=n, replace=False)] @ cell if n else None

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")

    plot_panel(
        ax1, cell, inv_cell, blocked_topo,
        f"Topological (probe={grid_topo.probe_radius} A)\nclosed pores in perfect crystal",
        t_pos=t_pos, t_edges=t_edges,
    )
    plot_panel(
        ax2, cell, inv_cell, blocked_steric,
        f"Steric (probe={grid_steric.probe_radius} A)\nconfined cages / narrow windows",
        t_pos=t_pos, t_edges=t_edges, accessible_cart=acc_sample,
    )
    fig.suptitle("FAU pore blocking + T-atom framework skeleton", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=160)
    plt.close(fig)
    print(f"[INFO] Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
