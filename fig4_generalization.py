#!/usr/bin/env python3
"""Figure 4: cross-framework generalization (zeolites, MOFs, COFs).

Caption summary: Accessibility-aware packing succeeds across six frameworks
and four guest species with zero violations; Packmol misplaces guests in closed
cages on cubic hosts.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from figstyle import use_nature_style, COLORS, panel_label, save_figure, DOUBLE_COL
from figutil import imshow_panel


def main():
    use_nature_style()
    d = np.load("runs/generalization.npz", allow_pickle=True)
    fws = [str(x) for x in d["frameworks"]]
    guests = [str(x) for x in d["guests"]]
    loading = d["loading"]
    viol = d["violations"]
    acc = d["accessible_cell_fraction"]
    pk_inacc = d["packmol_inaccessible"]

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.78))
    gs = GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32,
                  left=0.08, right=0.98, top=0.93, bottom=0.08,
                  height_ratios=[1.05, 0.95])

    heroes = [
        ("figures/panels/gen_mof5.png", "MOF-5 (MOF)"),
        ("figures/panels/gen_uio66.png", "UiO-66 (MOF)"),
        ("figures/panels/gen_cof5.png", "COF-5 (COF)"),
    ]
    for j, (png, title) in enumerate(heroes):
        ax = fig.add_subplot(gs[0, j])
        imshow_panel(ax, png)
        ax.set_title(title, fontsize=6.5, pad=3)
        panel_label(ax, chr(ord("a") + j), x=-0.04)

    # (d) loading heatmap
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(loading, aspect="auto", cmap="Blues")
    ax.set_xticks(range(len(guests)))
    ax.set_xticklabels(guests, rotation=45, ha="right")
    ax.set_yticks(range(len(fws)))
    ax.set_yticklabels(fws)
    for i in range(len(fws)):
        for j in range(len(guests)):
            txt = f"{loading[i, j]}"
            if viol[i, j] == 0:
                txt += "\n0v"
            else:
                txt += f"\n{viol[i, j]}v"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.5,
                    color="white" if loading[i, j] > loading.max() * 0.55 else "black")
    ax.set_xlabel("Guest")
    ax.set_ylabel("Framework")
    ax.set_title("Loading (molecules / cell)", fontsize=6.5)
    panel_label(ax, "d")

    # (e) accessible volume fraction
    ax = fig.add_subplot(gs[1, 1])
    colors = [COLORS["accessible"] if f == "zeolite" else
              ("#756bb1" if f == "MOF" else "#d95f0e")
              for f in d["families"]]
    ax.bar(range(len(fws)), 100 * acc, color=colors, width=0.7)
    ax.set_xticks(range(len(fws)))
    ax.set_xticklabels(fws, rotation=45, ha="right")
    ax.set_ylabel("Accessible volume (% of cell)")
    ax.set_title("Pore accessibility", fontsize=6.5)
    panel_label(ax, "e")

    # (f) Packmol closed-cage errors (cubic hosts only)
    ax = fig.add_subplot(gs[1, 2])
    cubic_idx = [i for i, fw in enumerate(fws) if fw in ("FAU", "MOF-5")]
    cubic_guests = guests
    x = np.arange(len(cubic_guests))
    w = 0.35
    for k, fi in enumerate(cubic_idx):
        vals = [0 if np.isnan(pk_inacc[fi, j]) else int(pk_inacc[fi, j])
                for j in range(len(cubic_guests))]
        offset = (k - 0.5) * w
        ax.bar(x + offset, vals, w, label=fws[fi],
               color=COLORS["ours"] if k == 0 else COLORS["packmol"])
    ax.axhline(0, color="k", lw=0.5)
    ax.plot(x, [0] * len(x), "D", color=COLORS["ours"], ms=4, label="Ours (0)")
    ax.set_xticks(x)
    ax.set_xticklabels(cubic_guests, rotation=45, ha="right")
    ax.set_ylabel("Guests in closed cages")
    ax.set_title("Packmol on cubic hosts", fontsize=6.5)
    ax.legend(fontsize=5, loc="upper right", handlelength=1.0)
    panel_label(ax, "f")

    pdf, png = save_figure(fig, "fig4_generalization")
    print("saved", pdf, png)


if __name__ == "__main__":
    main()
