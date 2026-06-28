#!/usr/bin/env python3
"""Shared helpers for composing figures (image insets, trimming)."""

import numpy as np
import matplotlib.image as mpimg


def imshow_panel(ax, png_path, title=None, title_color="black"):
    """Show a rendered PNG in an axis with whitespace trimmed; no ticks."""
    img = mpimg.imread(png_path)
    img = _trim_white(img)
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, color=title_color, pad=2)


def _trim_white(img, thresh=0.985):
    """Crop uniform white/transparent border around a rendered panel."""
    if img.ndim != 3:
        return img
    if img.shape[2] == 4:
        alpha = img[..., 3]
        mask = alpha > 0.01
    else:
        lum = img[..., :3].mean(axis=2)
        mask = lum < thresh
    if not mask.any():
        return img
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    pad = 6
    r0, r1 = max(0, rows[0] - pad), min(img.shape[0], rows[-1] + pad)
    c0, c1 = max(0, cols[0] - pad), min(img.shape[1], cols[-1] + pad)
    return img[r0:r1, c0:c1]
