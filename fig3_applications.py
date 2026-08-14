#!/usr/bin/env python3
"""Figure 3: electrolyte, CO2 adsorption, and separation vs Packmol.

Caption summary: Three application cases show zero physical violations for
accessibility-aware packing vs systematic closed-cage and clash errors in
Packmol.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from figstyle import use_nature_style, COLORS, panel_label, save_figure, DOUBLE_COL
from figutil import imshow_panel

CASES = ["electrolyte", "co2", "separation"]
TITLES = {"electrolyte": "Li-salt electrolyte (EC/DMC/LiPF$_6$)",
          "co2": "CO$_2$ adsorption",
          "separation": "CO$_2$/H$_2$O separation"}


def main():
    use_nature_style()
    data = {c: np.load(f"runs/app_{c}.npz", allow_pickle=True) for c in CASES}

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.95))
    gs = GridSpec(3, 3, figure=fig, hspace=0.34, wspace=0.32,
                  left=0.10, right=0.985, top=0.94, bottom=0.07,
                  height_ratios=[1.0, 1.0, 0.85])

    # ---- rows 0,1: renders (ours / packmol) ----
    for j, c in enumerate(CASES):
        ax = fig.add_subplot(gs[0, j])
        imshow_panel(ax, f"figures/panels/app_{c}_ours.png")
        ax.set_title(TITLES[c], fontsize=6.2, pad=3)
        if j == 0:
            ax.set_ylabel("Accessibility-aware\n(ours)", fontsize=6.5,
                          color=COLORS["ours"], labelpad=2)
            panel_label(ax, "a", x=-0.02)

        ax = fig.add_subplot(gs[1, j])
        imshow_panel(ax, f"figures/panels/app_{c}_packmol.png")
        if j == 0:
            ax.set_ylabel("Packmol\n(blind; errors in red)", fontsize=6.5,
                          color=COLORS["packmol"], labelpad=4)
            panel_label(ax, "b", x=-0.02)

    # ---- (c) inaccessible-cage errors per case ----
    ax = fig.add_subplot(gs[2, 0])
    x = np.arange(len(CASES))
    w = 0.36
    ours_inacc = [int(data[c]["ours_inacc"]) for c in CASES]
    pk_inacc = [int(data[c]["pk_inacc"]) for c in CASES]
    ax.bar(x - w / 2, ours_inacc, w, color=COLORS["ours"], label="Ours")
    ax.bar(x + w / 2, pk_inacc, w, color=COLORS["packmol"], label="Packmol")
    ax.set_xticks(x)
    ax.set_xticklabels(["elec.", "CO$_2$", "sep."])
    ax.set_ylabel("Guests in closed cages")
    ax.legend(loc="upper left", handlelength=1.2)
    panel_label(ax, "c")

    # ---- (d) all violation types (Packmol), stacked; ours overlaid ----
    ax = fig.add_subplot(gs[2, 1])
    inacc = np.array([int(data[c]["pk_inacc"]) for c in CASES])
    gf = np.array([int(data[c]["pk_gf"]) for c in CASES])
    gg = np.array([int(data[c]["pk_gg"]) for c in CASES])
    ax.bar(x, inacc, w + 0.1, color=COLORS["packmol"], label="Closed cage")
    ax.bar(x, gf, w + 0.1, bottom=inacc, color="#fc8d59", label="Framework clash")
    ax.bar(x, gg, w + 0.1, bottom=inacc + gf, color="#d9a300", label="Guest overlap")
    ours_tot = [int(data[c]["ours_inacc"]) + int(data[c]["ours_gf"]) +
                int(data[c]["ours_gg"]) for c in CASES]
    ax.plot(x, ours_tot, "D", color=COLORS["ours"], ms=4, label="Ours (total)")
    ax.set_xticks(x)
    ax.set_xticklabels(["elec.", "CO$_2$", "sep."])
    ax.set_ylabel("Physical violations (Packmol)")
    ax.legend(loc="upper left", handlelength=1.2)
    ax.set_title("All violation types", fontsize=6.5)
    panel_label(ax, "d")

    # ---- (e) guest-framework min-distance distribution ----
    ax = fig.add_subplot(gs[2, 2])
    ours_gf = np.concatenate([data[c]["ours_gf_min"] for c in CASES])
    pk_gf = np.concatenate([data[c]["pk_gf_min"] for c in CASES])
    bins = np.linspace(0, 6, 25)
    ax.hist(pk_gf, bins=bins, color=COLORS["packmol"], alpha=0.6,
            label="Packmol", density=True)
    ax.hist(ours_gf, bins=bins, color=COLORS["ours"], alpha=0.6,
            label="Ours", density=True)
    ax.axvline(2.2, color="k", ls=":", lw=0.8)
    ax.text(2.2, ax.get_ylim()[1] * 0.9, " cutoff", fontsize=5, va="top")
    ax.set_xlabel("Guest–framework min. distance (Å)")
    ax.set_ylabel("Density")
    ax.legend(loc="upper left", handlelength=1.2)
    panel_label(ax, "e")

    pdf, png = save_figure(fig, "fig3_applications")
    print("saved", pdf, png)


if __name__ == "__main__":
    main()
