#!/usr/bin/env python3
"""Figure 1: accessibility-aware packing workflow and FAU pore map.

Caption summary: Four-step workflow from host framework to validated host-guest
model; steric flood-fill marks inaccessible sodalite cages invisible to
geometry-only packers.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch

from figstyle import use_nature_style, COLORS, panel_label, save_figure, DOUBLE_COL
from figutil import imshow_panel
from pore_accessibility import PoreGrid

STERIC = "frameworks/zeolites/FAU.blocked.npz"

WORKFLOW = [
    ("figures/panels/fw_only.png", "1  Host framework\n+ guest library"),
    ("figures/panels/fw_steric_blocked.png", "2  Pore accessibility\n(flood-fill)"),
    ("figures/panels/packed_on.png", "3  Annealed MC\ninsertion"),
    ("figures/panels/packed_validated.png", "4  Validated\nhost–guest model"),
]


def voxel_slice(grid):
    """Categorical z-slice: 0 solid, 1 accessible, 2 blocked."""
    cat = np.zeros(grid.solid_mask.shape, dtype=int)
    cat[grid.accessible_mask] = 1
    cat[grid.blocked_mask] = 2
    # choose the slice with the most blocked voxels
    z = int(np.argmax(grid.blocked_mask.sum(axis=(0, 1))))
    return cat[:, :, z], z


def main():
    use_nature_style()
    grid = PoreGrid.load(STERIC)

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.68))
    gs = GridSpec(2, 4, figure=fig, hspace=0.12, wspace=0.28,
                  left=0.04, right=0.985, top=0.94, bottom=0.06,
                  height_ratios=[1.0, 1.15])

    # ---------- (a) workflow strip ----------
    axes_wf = []
    for j, (png, label) in enumerate(WORKFLOW):
        ax = fig.add_subplot(gs[0, j])
        imshow_panel(ax, png)
        ax.set_title(label, fontsize=6.0, pad=2)
        axes_wf.append(ax)
    panel_label(axes_wf[0], "a", x=-0.02, y=1.12)
    # arrows between workflow steps
    for j in range(len(axes_wf) - 1):
        a0 = axes_wf[j]
        a1 = axes_wf[j + 1]
        arr = FancyArrowPatch(
            (1.02, 0.5), (1.13, 0.5), transform=a0.transAxes,
            arrowstyle="-|>", mutation_scale=8, lw=1.0, color="0.3",
            clip_on=False,
        )
        fig.add_artist(arr)

    # ---------- (b) accessibility hero ----------
    axb = fig.add_subplot(gs[1, 0:2])
    imshow_panel(axb, "figures/panels/fw_steric_blocked.png")
    axb.set_title("Inaccessible sodalite cages in FAU", fontsize=6.5, pad=6)
    panel_label(axb, "b", x=-0.02, y=1.14)

    # ---------- (c) voxel cross-section ----------
    axc = fig.add_subplot(gs[1, 2])
    sl, z = voxel_slice(grid)
    cmap = ListedColormap([COLORS["solid"], COLORS["accessible"], COLORS["blocked"]])
    axc.imshow(sl.T, origin="lower", cmap=cmap, vmin=0, vmax=2, interpolation="nearest")
    axc.set_xticks([])
    axc.set_yticks([])
    axc.set_title("Voxel map (z-slice)", fontsize=6.5, pad=3)
    # legend
    from matplotlib.patches import Patch
    handles = [Patch(color=COLORS["solid"], label="Solid"),
               Patch(color=COLORS["accessible"], label="Accessible"),
               Patch(color=COLORS["blocked"], label="Blocked")]
    axc.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, fontsize=5, handlelength=1.0, columnspacing=0.8)
    panel_label(axc, "c")

    # ---------- (d) volume partition ----------
    axd = fig.add_subplot(gs[1, 3])
    s = grid.stats
    total = grid.grid_n ** 3
    vals = [100 * s["solid_voxels"] / total,
            100 * s["accessible_voxels"] / total,
            100 * s["blocked_voxels"] / total]
    cols = [COLORS["solid"], COLORS["accessible"], COLORS["blocked"]]
    bars = axd.bar(["Solid", "Access.", "Blocked"], vals, color=cols, width=0.7)
    for b, v in zip(bars, vals):
        axd.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.0f}",
                 ha="center", fontsize=5.5)
    axd.set_ylabel("Cell volume (%)")
    axd.set_title("FAU partition", fontsize=6.5, pad=3)
    axd.set_ylim(0, max(vals) * 1.18)
    panel_label(axd, "d")

    pdf, png = save_figure(fig, "fig1_overview")
    print("saved", pdf, png)


if __name__ == "__main__":
    main()
