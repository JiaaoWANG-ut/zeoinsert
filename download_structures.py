#!/usr/bin/env python3
"""Download zeolite / MOF / COF frameworks and small-molecule XYZ files."""

import math
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FRAMEWORKS = ROOT / "frameworks"
MOLECULES = ROOT / "molecules"

IZA_BASE = "https://www.iza-structure.org/IZA-SC/cif/{code}.cif"
WMD_BASE = (
    "https://raw.githubusercontent.com/WMD-group/Crystal_structures/master/{path}"
)
COF_BASE = (
    "https://raw.githubusercontent.com/danieleongari/CURATED-COFs/master/cifs/{fname}"
)
PUBCHEM_SDF = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF"
    "?record_type=3d"
)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  -> {dest.name}")
    req = urllib.request.Request(url, headers={"User-Agent": "zeoinsert/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()
    dest.write_bytes(data)


def sdf_to_xyz(sdf_text: str) -> str:
    lines = sdf_text.splitlines()
    counts = lines[3].split()
    n_atoms = int(counts[0])
    atoms = []
    for i in range(4, 4 + n_atoms):
        parts = lines[i].split()
        x, y, z = map(float, parts[:3])
        sym = re.sub(r"\d+", "", parts[3])
        atoms.append((sym, x, y, z))
    body = "\n".join(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}" for sym, x, y, z in atoms)
    return f"{n_atoms}\n\n{body}\n"


def write_xyz(path: Path, comment: str, atoms: list[tuple[str, tuple[float, float, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment]
    for sym, (x, y, z) in atoms:
        lines.append(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    path.write_text("\n".join(lines) + "\n")
    print(f"  -> {path.name}")


def fetch_pubchem_xyz(name: str, dest: Path) -> bool:
    url = PUBCHEM_SDF.format(name=urllib.parse.quote(name))
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "zeoinsert/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            sdf = resp.read().decode()
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(sdf_to_xyz(sdf))
        print(f"  -> {dest.name} (PubChem)")
        return True
    except urllib.error.HTTPError:
        return False


def build_manual_molecules() -> None:
    """Geometries from gas-phase / crystal ion-pair references (Å)."""
    pubchem_names = {
        "CO2.xyz": "carbon dioxide",
        "EC.xyz": "ethylene carbonate",
        "DMC.xyz": "dimethyl carbonate",
        "H2O.xyz": "water",
    }
    for fname, pname in pubchem_names.items():
        dest = MOLECULES / fname
        if not fetch_pubchem_xyz(pname, dest):
            print(f"  !! PubChem failed for {fname}, using fallback")

    # CO2 fallback (linear, C-O = 1.16 Å)
    if not (MOLECULES / "CO2.xyz").exists():
        write_xyz(
            MOLECULES / "CO2.xyz",
            "CO2 linear",
            [("C", (0.0, 0.0, 0.0)), ("O", (-1.16, 0.0, 0.0)), ("O", (1.16, 0.0, 0.0))],
        )

    # Li+ (single ion, for GCMC / insertion tests)
    write_xyz(MOLECULES / "Li_ion.xyz", "Li+ single ion", [("Li", (0.0, 0.0, 0.0))])

    # LiCl ion pair (~2.02 Å)
    write_xyz(
        MOLECULES / "LiCl.xyz",
        "LiCl ion pair",
        [("Li", (0.0, 0.0, 0.0)), ("Cl", (2.02, 0.0, 0.0))],
    )

    # LiPF6: octahedral PF6- with Li+ at ~2.5 Å along z
    pf_bond = 1.58
    li_p = 2.50
    atoms_pf6 = [("P", (0.0, 0.0, 0.0))]
    for i in range(6):
        theta = math.acos(1.0 - 2.0 * (i + 0.5) / 6.0)
        phi = math.pi * (1.0 + 5.0 ** 0.5) * i
        atoms_pf6.append(
            (
                "F",
                (
                    pf_bond * math.sin(theta) * math.cos(phi),
                    pf_bond * math.sin(theta) * math.sin(phi),
                    pf_bond * math.cos(theta),
                ),
            )
        )
    atoms_lipf6 = [("Li", (0.0, 0.0, -li_p))] + atoms_pf6
    write_xyz(MOLECULES / "LiPF6.xyz", "LiPF6 ion pair (approx.)", atoms_lipf6)

    # LiBF4: tetrahedral BF4- with Li+ at ~2.0 Å
    bf_bond = 1.38
    li_b = 2.00
    tet_dirs = [
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
    ]
    atoms_bf4 = [("B", (0.0, 0.0, 0.0))]
    for dx, dy, dz in tet_dirs:
        norm = (dx * dx + dy * dy + dz * dz) ** 0.5
        atoms_bf4.append(("F", (bf_bond * dx / norm, bf_bond * dy / norm, bf_bond * dz / norm)))
    atoms_libf4 = [("Li", (0.0, 0.0, -li_b))] + atoms_bf4
    write_xyz(MOLECULES / "LiBF4.xyz", "LiBF4 ion pair (approx.)", atoms_libf4)

    # LiTFSI simplified geometry (N at center, two -SO2CF3 arms)
    # Based on typical gas-phase conformer (approximate)
    write_xyz(
        MOLECULES / "LiTFSI.xyz",
        "LiTFSI approximate conformer",
        [
            ("Li", (0.00, 0.00, -2.80)),
            ("N", (0.00, 0.00, 0.00)),
            ("S", (-1.55, 0.00, 1.10)),
            ("S", (1.55, 0.00, 1.10)),
            ("O", (-2.35, 1.10, 0.55)),
            ("O", (-2.35, -1.10, 0.55)),
            ("O", (2.35, 1.10, 0.55)),
            ("O", (2.35, -1.10, 0.55)),
            ("C", (-1.55, 0.00, 2.65)),
            ("C", (1.55, 0.00, 2.65)),
            ("F", (-1.55, 0.95, 3.45)),
            ("F", (-1.55, -0.95, 3.45)),
            ("F", (-2.50, 0.00, 2.65)),
            ("F", (1.55, 0.95, 3.45)),
            ("F", (1.55, -0.95, 3.45)),
            ("F", (2.50, 0.00, 2.65)),
        ],
    )


def main() -> None:
    print("=== Frameworks (CIF) ===\n")

    print("[Zeolites]")
    zeolites = {
        "FAU.cif": "FAU",           # cubic (reference FAU)
        "HEU.cif": "HEU",           # monoclinic, beta=116°
        "ERI.cif": "ERI",           # hexagonal, gamma=120°
        "LTL.cif": "LTL",           # hexagonal, gamma=120°
        "OFF.cif": "OFF",           # hexagonal, gamma=120°
        "MAZ.cif": "MAZ",           # hexagonal, gamma=120°
    }
    for fname, code in zeolites.items():
        download(IZA_BASE.format(code=code), FRAMEWORKS / "zeolites" / fname)

    print("\n[MOFs]")
    mofs = {
        "MOF-5.cif": "MOFs/MOF-5/MOF5.cif",
        "UiO-66.cif": "MOFs/UiO/UiO-66.cif",
        "HKUST-1.cif": "MOFs/HKUST-1/HKUST1.cif",
        "MIL-125.cif": "MOFs/MIL-125/2014_PBEsol/MIL125.cif",
    }
    for fname, path in mofs.items():
        download(WMD_BASE.format(path=path), FRAMEWORKS / "mofs" / fname)

    print("\n[COFs]")
    cofs = {
        "COF-1.cif": "05000N2.cif",
        "COF-5.cif": "05001N2.cif",
        "COF-102.cif": "07010N3.cif",
        "COF-202.cif": "08000N3.cif",
    }
    for fname, src in cofs.items():
        download(COF_BASE.format(fname=src), FRAMEWORKS / "cofs" / fname)

    print("\n=== Small molecules (XYZ) ===\n")
    build_manual_molecules()

    print("\nDone.")
    print(f"  Frameworks: {FRAMEWORKS}")
    print(f"  Molecules:  {MOLECULES}")


if __name__ == "__main__":
    main()
