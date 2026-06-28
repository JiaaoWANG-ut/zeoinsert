#!/usr/bin/env python3
"""Figure 2: probe physics, ablation vs blind packing, and MC convergence.

Caption summary: Guest-size-dependent pore closure, zero misplacement with
accessibility-aware packing, and simulated-annealing convergence.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from figstyle import use_nature_style, COLORS, panel_label, save_figure, DOUBLE_COL
from figutil import imshow_panel


def main():
    use_nature_style()
    sweep = np.load("runs/probe_sweep.npz", allow_pickle=True)
    abl = np.load("runs/ablation.npz")
    conv = np.load("runs/convergence.npz")

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.62))
    gs = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.40,
                  left=0.07, right=0.985, top=0.93, bottom=0.10)

    # ---------- (a) probe-radius sweep ----------
    ax = fig.add_subplot(gs[0, 0])
    radii = sweep["radii"]
    fws = [str(x) for x in sweep["frameworks"]]
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, len(fws)))
    for c, f in zip(cmap, fws):
        ax.plot(radii, 100 * sweep[f"{f}_acc_cell"], "-", color=c, label=f, lw=1.1)
    ax.set_xlabel("Probe radius (Å)")
    ax.set_ylabel("Accessible volume (% of cell)")
    ax.set_xlim(0.4, float(radii.max()) + 0.15)
    ax.legend(ncol=1, loc="upper right", handlelength=1.0, labelspacing=0.2,
              fontsize=5.5)
    # guest markers
    gnames = [str(x) for x in sweep["guest_names"]]
    grad = sweep["guest_radii"]
    ymax = ax.get_ylim()[1]
    for nm, rr in zip(gnames, grad):
        ax.axvline(rr, color="0.55", ls=":", lw=0.6, zorder=0)
        ax.text(rr + 0.04, ymax * 0.97, nm, rotation=90, va="top", ha="left",
                fontsize=5, color="0.35", clip_on=True)
    panel_label(ax, "a")

    # ---------- (b) blocked-fraction step (physics) ----------
    ax = fig.add_subplot(gs[0, 1])
    for c, f in zip(cmap, fws):
        ax.plot(radii, 100 * sweep[f"{f}_blocked_void"], "-", color=c, lw=1.1, label=f)
    ax.set_xlabel("Probe radius (Å)")
    ax.set_ylabel("Inaccessible void (%)")
    ax.set_title("Cage closure vs guest size", fontsize=6.5)
    panel_label(ax, "b")

    # ---------- (c) ablation: misplacement vs loading ----------
    ax = fig.add_subplot(gs[0, 2])
    load = abl["loadings"]
    mon = 100 * abl["mis_on"]
    moff = 100 * abl["mis_off"]
    ax.plot(load, mon.mean(1), "o-", color=COLORS["ours"], label="Accessibility-aware (ours)")
    ax.fill_between(load, mon.mean(1) - mon.std(1), mon.mean(1) + mon.std(1),
                    color=COLORS["ours"], alpha=0.2)
    ax.plot(load, moff.mean(1), "s--", color=COLORS["packmol"], label="Accessibility-blind")
    ax.fill_between(load, moff.mean(1) - moff.std(1), moff.mean(1) + moff.std(1),
                    color=COLORS["packmol"], alpha=0.2)
    ax.set_xlabel("Loading (molecules / cell)")
    ax.set_ylabel("Guests in closed cages (%)")
    ax.legend(loc="upper left", handlelength=1.5)
    ax.set_ylim(bottom=-0.5)
    panel_label(ax, "c")

    # ---------- (d,e) structure contrast ON vs OFF ----------
    axd = fig.add_subplot(gs[1, 0])
    imshow_panel(axd, "figures/panels/packed_on.png",
                 title="Accessibility-aware (ours)", title_color=COLORS["ours"])
    panel_label(axd, "d", x=-0.05)
    axe = fig.add_subplot(gs[1, 1])
    imshow_panel(axe, "figures/panels/packed_off.png",
                 title="Accessibility-blind (misplaced in red)", title_color=COLORS["packmol"])
    panel_label(axe, "e", x=-0.05)

    # ---------- (f) convergence ----------
    ax = fig.add_subplot(gs[1, 2])
    steps = conv["steps"]
    E = conv["energy"]
    T = conv["temperature"]
    meanE = np.nanmean(E, axis=0)
    valid = ~np.isnan(meanE)
    ax.plot(steps[valid], meanE[valid], color=COLORS["accessible"], lw=1.2,
            label="Overlap energy")
    for i in range(E.shape[0]):
        ax.plot(steps, E[i], color=COLORS["accessible"], lw=0.4, alpha=0.25)
    ax.set_xlabel("MC step")
    ax.set_ylabel("Overlap energy (a.u.)")
    ax.set_yscale("symlog", linthresh=0.1)
    axt = ax.twinx()
    axt.plot(steps, np.nanmean(T, axis=0), color="0.55", lw=0.9, ls="--")
    axt.set_ylabel("Temperature", color="0.45")
    axt.tick_params(axis="y", colors="0.45")
    axt.spines["right"].set_visible(True)
    axt.spines["top"].set_visible(False)
    ax.set_title("Simulated-annealing convergence", fontsize=6.5)
    panel_label(ax, "f")

    pdf, png = save_figure(fig, "fig2_probe_ablation")
    print("saved", pdf, png)


if __name__ == "__main__":
    main()
