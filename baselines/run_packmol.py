#!/usr/bin/env python3
"""Run Packmol (real third-party packer) as an accessibility-blind baseline.

Packmol enforces a global minimum inter-atom distance (tolerance) against the
fixed framework and between guests, but has no concept of pore accessibility,
so it places guests inside sterically closed cages. We feed it the same cell,
loadings, and cutoff as our method for a fair comparison.

Only orthogonal cells are supported (FAU, MOF-5, UiO-66 are cubic).
"""

from __future__ import annotations

import os
import subprocess
import tempfile

import numpy as np
from ovito.io import import_file

PACKMOL_BIN = os.environ.get("PACKMOL_BIN", "packmol")


def _get_symbols(data):
    if "Element" in data.particles:
        arr = data.particles["Element"].array
        return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]
    type_map = {t.id: t.name for t in data.particles.particle_types.types}
    return [type_map[i] for i in data.particles["Particle Type"].array]


def _read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.array(coords)


def _write_xyz(path, syms, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, c in zip(syms, coords):
            f.write(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")


def packmol_pack(framework_file, species, cell, tolerance=2.2, inset=2.5, seed=12345):
    """Pack guests into the framework cell with Packmol.

    species : list of dicts {name, file, count}
    cell    : 3x3 (must be ~orthogonal)
    Returns dict with combined positions/symbols, n_framework, per-molecule
    positions list, centers_frac, guest_atom_counts, guest_type_ids.
    """
    off = np.abs(cell - np.diag(np.diag(cell)))
    if off.max() > 1e-3:
        raise ValueError("Packmol baseline supports orthogonal cells only.")
    a, b, c = float(cell[0, 0]), float(cell[1, 1]), float(cell[2, 2])

    fw = import_file(framework_file).compute()
    fw_pos = np.array(fw.particles.positions)
    fw_sym = _get_symbols(fw)
    inv_cell = np.linalg.inv(cell)

    workdir = tempfile.mkdtemp(prefix="packmol_", dir="runs")
    fw_xyz = os.path.join(workdir, "framework.xyz")
    _write_xyz(fw_xyz, fw_sym, fw_pos, "framework")

    out_xyz = os.path.join(workdir, "packed.xyz")

    inp = [
        f"tolerance {tolerance}",
        "filetype xyz",
        f"seed {seed}",
        f"output {out_xyz}",
        "",
        f"structure {fw_xyz}",
        "  number 1",
        "  fixed 0. 0. 0. 0. 0. 0.",
        "end structure",
        "",
    ]
    guest_type_ids, guest_atom_counts, species_names = [], [], []
    for tid, sp in enumerate(species):
        count = int(sp["count"])
        if count <= 0:
            continue
        gsym, gcoord = _read_xyz(sp["file"])
        gfile = os.path.join(workdir, f"guest_{tid}.xyz")
        _write_xyz(gfile, gsym, gcoord, sp.get("name", ""))
        inp += [
            f"structure {gfile}",
            f"  number {count}",
            f"  inside box {inset:.3f} {inset:.3f} {inset:.3f} "
            f"{a - inset:.3f} {b - inset:.3f} {c - inset:.3f}",
            "end structure",
            "",
        ]
        guest_type_ids += [tid] * count
        guest_atom_counts += [len(gsym)] * count
        species_names.append(sp.get("name", f"sp{tid}"))

    inp_path = os.path.join(workdir, "pack.inp")
    with open(inp_path, "w") as f:
        f.write("\n".join(inp))

    with open(inp_path) as fin:
        proc = subprocess.run([PACKMOL_BIN], stdin=fin, capture_output=True,
                              text=True, timeout=600)
    if not os.path.exists(out_xyz):
        raise RuntimeError(f"Packmol failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-1000:]}")

    all_sym, all_pos = _read_xyz(out_xyz)
    n_fw = len(fw_sym)
    guest_pos = all_pos[n_fw:]

    mol_positions, centers, start = [], [], 0
    for c_ in guest_atom_counts:
        m = guest_pos[start:start + c_]
        mol_positions.append(m)
        centers.append(m.mean(axis=0))
        start += c_
    centers = np.array(centers) if centers else np.zeros((0, 3))
    centers_frac = (centers @ inv_cell) % 1.0 if len(centers) else centers

    return {
        "positions": all_pos,
        "symbols": all_sym,
        "cell": cell,
        "n_framework": n_fw,
        "mol_positions": mol_positions,
        "centers_frac": centers_frac,
        "guest_atom_counts": guest_atom_counts,
        "guest_type_ids": np.asarray(guest_type_ids, dtype=int),
        "species_names": species_names,
        "stdout": proc.stdout,
    }
