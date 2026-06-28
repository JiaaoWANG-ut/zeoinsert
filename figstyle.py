#!/usr/bin/env python3
"""Nature-style matplotlib configuration, palette, and panel helpers."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Nature column widths (mm -> inch)
MM = 1.0 / 25.4
SINGLE_COL = 89 * MM
DOUBLE_COL = 183 * MM

# consistent palette
COLORS = {
    "framework": "#9aa0a6",
    "accessible": "#2c7fb8",
    "blocked": "#d7301f",
    "solid": "#cfd4d9",
    "ours": "#1a9850",       # our method (good)
    "packmol": "#d73027",    # external baseline (errors)
    "random": "#8073ac",     # naive baseline
    "co2": "#e6550d",
    "h2o": "#3182bd",
    "li": "#9467bd",
    "ec": "#31a354",
    "dmc": "#756bb1",
    "violation": "#d7301f",
}

SPECIES_COLORS = {
    "CO2": "#e6550d", "H2O": "#3182bd", "EC": "#31a354",
    "DMC": "#756bb1", "LiPF6": "#9467bd", "Li": "#9467bd",
}


def use_nature_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,       # editable text in PDF
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def panel_label(ax, letter, x=-0.18, y=1.05, fontsize=9):
    """Bold lowercase panel letter, Nature style."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="right")


def save_figure(fig, stem, outdir="figures"):
    import os
    os.makedirs(outdir, exist_ok=True)
    pdf = os.path.join(outdir, f"{stem}.pdf")
    png = os.path.join(outdir, f"{stem}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    plt.close(fig)
    return pdf, png
