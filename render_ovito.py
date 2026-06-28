#!/usr/bin/env python3
"""High-quality OVITO Tachyon rendering for host-guest structures.

Framework is rendered as a muted grey skeleton (small spheres + bonds); guests
are vivid CPK spheres. A per-molecule `bad_mask` highlights violating guests in
red so that comparison panels (e.g. Packmol vs ours) read at a glance.
"""

from __future__ import annotations

import os

import numpy as np
from ovito.data import DataCollection, Particles, SimulationCell, Bonds
from ovito.pipeline import Pipeline, StaticSource
from ovito.modifiers import CreateBondsModifier
from ovito.vis import Viewport, TachyonRenderer, BondsVis  # noqa: F401

# CPK-ish colors for guest atoms
ELEM_COLOR = {
    "H": (0.95, 0.95, 0.95), "C": (0.20, 0.20, 0.20), "N": (0.19, 0.31, 0.97),
    "O": (1.00, 0.05, 0.05), "F": (0.56, 0.88, 0.31), "P": (1.00, 0.50, 0.00),
    "S": (1.00, 1.00, 0.19), "Cl": (0.12, 0.94, 0.12), "Li": (0.80, 0.50, 1.00),
    "B": (1.00, 0.71, 0.71), "Si": (0.94, 0.78, 0.63), "Al": (0.75, 0.65, 0.65),
}
ELEM_RADIUS = {
    "H": 0.30, "C": 0.55, "N": 0.55, "O": 0.55, "F": 0.50, "P": 0.72,
    "S": 0.66, "Cl": 0.66, "Li": 0.58, "B": 0.55, "Si": 0.55, "Al": 0.55,
}
FRAMEWORK_COLOR = (0.72, 0.74, 0.78)
FRAMEWORK_RADIUS = 0.40
HIGHLIGHT_COLOR = (0.84, 0.19, 0.12)
BLOCKED_COLOR = (0.84, 0.19, 0.12)


def _elem_color(sym):
    return ELEM_COLOR.get(sym, (0.5, 0.5, 0.5))


def _elem_radius(sym):
    return ELEM_RADIUS.get(sym, 0.55)


def render_structure(
    positions,
    symbols,
    cell,
    n_framework,
    out_png,
    guest_atom_counts=None,
    bad_mask=None,
    blocked_cart=None,
    blocked_radius=0.6,
    blocked_radii=None,
    blocked_max=1500,
    blocked_transparency=0.55,
    size=(900, 900),
    camera_dir=(-1.0, -1.2, -0.8),
    fov=None,
    framework_bonds=True,
    bond_cutoff=2.0,
    background=(1, 1, 1),
    alpha=True,
    ambient_occlusion=True,
):
    """Render a host-guest structure to PNG.

    bad_mask : optional (n_guest_molecules,) bool; flagged molecules drawn red.
    blocked_cart : optional (M,3) cartesian centers of blocked cages (translucent).
    """
    positions = np.asarray(positions, dtype=float)
    n_total = len(positions)
    n_guest = n_total - n_framework

    colors = np.zeros((n_total, 3))
    radii = np.zeros(n_total)
    transp = np.zeros(n_total)

    # framework
    colors[:n_framework] = FRAMEWORK_COLOR
    radii[:n_framework] = FRAMEWORK_RADIUS

    # guests
    guest_syms = symbols[n_framework:]
    for i, s in enumerate(guest_syms):
        colors[n_framework + i] = _elem_color(s)
        radii[n_framework + i] = _elem_radius(s)

    # highlight violating molecules
    if bad_mask is not None and guest_atom_counts is not None and n_guest > 0:
        start = n_framework
        for m, c in enumerate(guest_atom_counts):
            if m < len(bad_mask) and bad_mask[m]:
                colors[start:start + c] = HIGHLIGHT_COLOR
                radii[start:start + c] *= 1.15
            start += c

    data = DataCollection()
    sc = SimulationCell(pbc=(True, True, True))
    sc[...] = np.column_stack([cell.T, np.zeros(3)])
    sc.vis.enabled = True
    sc.vis.line_width = max(cell.diagonal()) * 0.004
    sc.vis.rendering_color = (0.0, 0.0, 0.0)
    data.objects.append(sc)

    parts = Particles()
    parts.create_property("Position", data=positions)
    parts.create_property("Color", data=colors)
    parts.create_property("Radius", data=radii)
    if transp.any():
        parts.create_property("Transparency", data=transp)
    parts.vis.enabled = True
    data.objects.append(parts)

    pipe = Pipeline(source=StaticSource(data=data))

    if framework_bonds:
        bonds_mod = CreateBondsModifier(cutoff=bond_cutoff)
        bonds_mod.vis.enabled = True
        bonds_mod.vis.width = 0.18
        bonds_mod.vis.flat_shading = False
        bonds_mod.vis.coloring_mode = BondsVis.ColoringMode.Uniform
        bonds_mod.vis.color = FRAMEWORK_COLOR
        pipe.modifiers.append(bonds_mod)

    pipe.add_to_scene()

    # optional translucent blocked-cage markers (separate pipeline)
    blocked_pipe = None
    if blocked_cart is not None and len(blocked_cart):
        bc = np.asarray(blocked_cart, dtype=float)
        if blocked_radii is not None:
            br = np.asarray(blocked_radii, dtype=float)
        else:
            if len(bc) > blocked_max:
                sel = np.random.default_rng(0).choice(len(bc), blocked_max, replace=False)
                bc = bc[sel]
            br = np.full(len(bc), blocked_radius)
        bdata = DataCollection()
        bparts = Particles()
        bparts.create_property("Position", data=bc)
        bparts.create_property("Color", data=np.tile(BLOCKED_COLOR, (len(bc), 1)))
        bparts.create_property("Radius", data=br)
        bparts.create_property("Transparency", data=np.full(len(bc), blocked_transparency))
        bparts.vis.enabled = True
        bdata.objects.append(bparts)
        blocked_pipe = Pipeline(source=StaticSource(data=bdata))
        blocked_pipe.add_to_scene()

    vp = Viewport(type=Viewport.Type.Perspective, camera_dir=camera_dir)
    if fov is not None:
        vp.fov = fov
    vp.zoom_all(size=size)

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    renderer = TachyonRenderer(shadows=False, ambient_occlusion=ambient_occlusion)
    vp.render_image(filename=out_png, size=size, renderer=renderer,
                    background=background, alpha=alpha)

    pipe.remove_from_scene()
    if blocked_pipe is not None:
        blocked_pipe.remove_from_scene()
    return out_png
