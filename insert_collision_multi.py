#!/usr/bin/env python3
"""Monte Carlo placement of multiple molecular species into a periodic framework."""

import numpy as np
from ovito.io import import_file

from pore_accessibility import PoreGrid

# ========== 用户参数 ==========
FRAMEWORK_FILE = "frameworks/zeolites/FAU.cif"
OUTPUT_FILE = "FAU_multi.cif"

N_TOTAL = 12

SPECIES = [
    {"name": "CO2", "file": "molecules/CO2.xyz", "ratio": 2},
    {"name": "EC", "file": "molecules/EC.xyz", "ratio": 1},
]

CUTOFF_FW = 2.2
CUTOFF_MM = 2.2

# ========== 封闭孔 / 可达性 ==========
USE_PORE_BLOCKING = True
BLOCKED_FILE = "frameworks/zeolites/FAU.blocked.npz"  # None = 运行时重新探测
GRID_N = 64
PROBE_RADIUS = 1.8          # Å，探针半径，建议 ≈ 最大插入物种半径
MANUAL_BLOCKED_BOXES = []     # 额外手动 block，分数坐标盒子
BLOCKED_PENALTY = 100.0       # 分子中心落入封闭孔时的罚分

MAX_ITERS = 30000
PRINT_EVERY = 200

STEP_TRANSL = 0.8
STEP_ROT_DEG = 20.0

T0 = 1.0
T_MIN = 0.02

P_SMALL_MOVE = 0.70
P_BIG_JUMP = 0.10
P_BIG_ROT = 0.10

E_TOL = 1e-6


# ========== 工具函数 ==========
def random_unit_vector():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)


def random_small_rotation(max_angle_deg):
    theta = np.deg2rad(max_angle_deg) * (2 * np.random.rand() - 1.0)
    axis = random_unit_vector()
    kx, ky, kz = axis
    k = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    i = np.eye(3)
    return i + np.sin(theta) * k + (1 - np.cos(theta)) * (k @ k)


def random_rotation():
    u1, u2, u3 = np.random.rand(3)
    q = np.array([
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3),
    ])
    w, x, y, z = q
    return np.array([
        [1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * z * w, 2 * x * z + 2 * y * w],
        [2 * x * y + 2 * z * w, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * x * w],
        [2 * x * z - 2 * y * w, 2 * y * z + 2 * x * w, 1 - 2 * x * x - 2 * y * y],
    ])


def get_symbols(data):
    if "Element" in data.particles:
        arr = data.particles["Element"].array
        return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]
    if "Particle Type" in data.particles:
        type_ids = data.particles["Particle Type"].array
        types = data.particles.particle_types.types
        type_map = {t.id: t.name for t in types}
        return [type_map[i] for i in type_ids]
    raise RuntimeError("No Element info found in data.")


def minimum_image(vecs, cell, inv_cell):
    frac = vecs @ inv_cell
    frac -= np.round(frac)
    return frac @ cell


def min_distance(a, b, cell, inv_cell):
    diff = a[:, None, :] - b[None, :, :]
    diff = minimum_image(diff, cell, inv_cell)
    return float(np.min(np.linalg.norm(diff, axis=2)))


def wrap_to_cell(pos, cell, inv_cell):
    frac = pos @ inv_cell
    frac -= np.floor(frac)
    return frac @ cell


def distribute_by_ratio(ratios, n_total):
    ratios = np.asarray(ratios, dtype=float)
    if np.any(ratios <= 0):
        raise ValueError("All ratios must be positive.")
    ideal = ratios / ratios.sum() * n_total
    counts = np.floor(ideal).astype(int)
    remainder = n_total - counts.sum()
    if remainder:
        order = np.argsort(-(ideal - counts))
        for k in range(remainder):
            counts[order[k % len(counts)]] += 1
    return counts


def build_species_list(species_cfg, n_total, default_cutoff_fw, default_cutoff_mm):
    has_count = all("count" in sp for sp in species_cfg)
    has_ratio = all("ratio" in sp for sp in species_cfg)

    if has_count:
        counts = [int(sp["count"]) for sp in species_cfg]
    elif has_ratio:
        if n_total is None or n_total <= 0:
            raise ValueError("N_TOTAL must be positive when using ratio mode.")
        counts = distribute_by_ratio([sp["ratio"] for sp in species_cfg], n_total).tolist()
    else:
        raise ValueError("Each species must define either 'ratio' or 'count' (all entries same mode).")

    if sum(counts) == 0:
        raise ValueError("Total molecule count is zero.")

    loaded = []
    mol_type_ids = []

    for idx, (sp, count) in enumerate(zip(species_cfg, counts)):
        if count <= 0:
            continue
        path = sp["file"]
        name = sp.get("name", path)
        print(f"[INFO] Loading species {idx} ({name}): {path}")
        mol = import_file(path).compute()
        pos = np.array(mol.particles.positions)
        pos -= pos.mean(axis=0)
        sym = get_symbols(mol)
        loaded.append({
            "name": name,
            "file": path,
            "pos": pos,
            "sym": sym,
            "cutoff_fw": sp.get("cutoff_fw", default_cutoff_fw),
            "cutoff_mm": sp.get("cutoff_mm", default_cutoff_mm),
            "count": count,
        })
        print(f"       atoms={len(pos)}, count={count}")
        mol_type_ids.extend([len(loaded) - 1] * count)

    return loaded, np.asarray(mol_type_ids, dtype=int)


def pair_cutoff_mm(species, type_i, type_j):
    return max(species[type_i]["cutoff_mm"], species[type_j]["cutoff_mm"])


def load_pore_grid(pos_fw, cell, inv_cell):
    if not USE_PORE_BLOCKING:
        return None
    if BLOCKED_FILE:
        print(f"[INFO] Loading blocked pore map: {BLOCKED_FILE}")
        grid = PoreGrid.load(BLOCKED_FILE)
    else:
        print("[INFO] Probing pore accessibility (no BLOCKED_FILE set)...")
        grid = PoreGrid.build(
            pos_fw,
            cell,
            inv_cell,
            grid_n=GRID_N,
            probe_radius=PROBE_RADIUS,
            manual_boxes=MANUAL_BLOCKED_BOXES,
            framework_file=FRAMEWORK_FILE,
        )
    grid.print_report()
    return grid


def random_center(pore_grid):
    if pore_grid is None:
        return np.random.rand(3)
    return pore_grid.random_accessible_center()


# ========== 读取框架 ==========
print(f"[INFO] Loading framework: {FRAMEWORK_FILE}")
fw = import_file(FRAMEWORK_FILE).compute()
pos_fw = np.array(fw.particles.positions)
sym_fw = get_symbols(fw)
cell = np.array(fw.cell.matrix)[:3, :3]
inv_cell = np.linalg.inv(cell)
print(f"[INFO] Framework atoms: {len(pos_fw)}")

pore_grid = load_pore_grid(pos_fw, cell, inv_cell)

species, mol_type_ids = build_species_list(SPECIES, N_TOTAL, CUTOFF_FW, CUTOFF_MM)
n_molecules = len(mol_type_ids)
print(f"\n[INFO] Total molecules to insert: {n_molecules}")
for sp in species:
    print(f"       {sp['name']}: {sp['count']}")

centers_frac = np.array([random_center(pore_grid) for _ in range(n_molecules)])
rot_mats = np.array([random_rotation() for _ in range(n_molecules)])


def world_positions(centers_frac_local, rot_mats_local):
    positions = []
    for i in range(n_molecules):
        sp = species[mol_type_ids[i]]
        center_cart = centers_frac_local[i] @ cell
        positions.append(sp["pos"] @ rot_mats_local[i].T + center_cart)
    return positions


def compute_conflicts_and_energy(world_list, centers_frac_local):
    conflicts = np.zeros(n_molecules)
    e_tot = 0.0

    if pore_grid is not None:
        for i in range(n_molecules):
            if pore_grid.is_center_blocked(centers_frac_local[i]):
                e_tot += BLOCKED_PENALTY
                conflicts[i] += BLOCKED_PENALTY

    for i in range(n_molecules):
        ti = mol_type_ids[i]
        d = min_distance(world_list[i], pos_fw, cell, inv_cell)
        cutoff = species[ti]["cutoff_fw"]
        if d < cutoff:
            e = (cutoff - d) ** 2
            e_tot += e
            conflicts[i] += e

    for i in range(n_molecules):
        for j in range(i + 1, n_molecules):
            d = min_distance(world_list[i], world_list[j], cell, inv_cell)
            cutoff = pair_cutoff_mm(species, mol_type_ids[i], mol_type_ids[j])
            if d < cutoff:
                e = (cutoff - d) ** 2
                e_tot += e
                conflicts[i] += e
                conflicts[j] += e

    return conflicts, e_tot


# ========== 模拟退火 MC ==========
print("\n[INFO] Starting global MC with simulated annealing (multi-species)...")

world_list = world_positions(centers_frac, rot_mats)
conflicts, e = compute_conflicts_and_energy(world_list, centers_frac)
print(f"[INFO] Initial E = {e:.4f}")

for it in range(1, MAX_ITERS + 1):
    t = max(T_MIN, T0 * (1.0 - it / MAX_ITERS))

    if e < E_TOL:
        print(f"\n✅ Converged at iter {it}, E={e:.6e}")
        break

    worst = int(np.argmax(conflicts))
    frac_old = centers_frac[worst].copy()
    r_old = rot_mats[worst].copy()

    r = np.random.rand()
    if r < P_SMALL_MOVE:
        d_cart = STEP_TRANSL * (2 * np.random.rand() - 1) * random_unit_vector()
        d_frac = d_cart @ inv_cell
        centers_frac[worst] = (frac_old + d_frac) % 1.0
        rot_mats[worst] = random_small_rotation(STEP_ROT_DEG) @ r_old
    elif r < P_SMALL_MOVE + P_BIG_JUMP:
        centers_frac[worst] = random_center(pore_grid)
        rot_mats[worst] = random_rotation()
    else:
        centers_frac[worst] = frac_old
        rot_mats[worst] = random_small_rotation(60.0) @ r_old

    world_new = world_positions(centers_frac, rot_mats)
    conflicts_new, e_new = compute_conflicts_and_energy(world_new, centers_frac)
    de = e_new - e

    if de <= 0.0 or np.random.rand() < np.exp(-de / t):
        world_list = world_new
        conflicts = conflicts_new
        e = e_new
    else:
        centers_frac[worst] = frac_old
        rot_mats[worst] = r_old

    if it % PRINT_EVERY == 0:
        sp_name = species[mol_type_ids[worst]]["name"]
        blocked_flag = (
            pore_grid.is_center_blocked(centers_frac[worst]) if pore_grid else False
        )
        print(
            f"[ITER {it}] E={e:.4f}, T={t:.3f}, "
            f"worst={worst}({sp_name}), blocked={blocked_flag}, "
            f"conflict={conflicts[worst]:.4f}"
        )
else:
    print(f"\n⚠ Reached MAX_ITERS={MAX_ITERS}, E={e:.4f} (may still have mild overlaps).")

if pore_grid is not None:
    n_blocked = sum(pore_grid.is_center_blocked(centers_frac[i]) for i in range(n_molecules))
    print(f"[INFO] Molecules in blocked pores after MC: {n_blocked}/{n_molecules}")

# ========== 拼接最终结构 ==========
all_positions = pos_fw.copy()
all_symbols = list(sym_fw)

for i in range(n_molecules):
    sp = species[mol_type_ids[i]]
    all_positions = np.vstack((all_positions, world_list[i]))
    all_symbols.extend(sp["sym"])

all_positions = wrap_to_cell(all_positions, cell, inv_cell)

# ========== 写 CIF ==========
print("\n[INFO] Writing final CIF:", OUTPUT_FILE)

a = np.linalg.norm(cell[0])
b = np.linalg.norm(cell[1])
c = np.linalg.norm(cell[2])
alpha = np.degrees(np.arccos(np.dot(cell[1], cell[2]) / (b * c)))
beta = np.degrees(np.arccos(np.dot(cell[0], cell[2]) / (a * c)))
gamma = np.degrees(np.arccos(np.dot(cell[0], cell[1]) / (a * b)))

with open(OUTPUT_FILE, "w") as f:
    f.write("data_mc_multi_sa\n")
    f.write("_symmetry_space_group_name_H-M 'P1'\n")
    f.write("_symmetry_Int_Tables_number 1\n\n")
    f.write(f"_cell_length_a {a:.6f}\n")
    f.write(f"_cell_length_b {b:.6f}\n")
    f.write(f"_cell_length_c {c:.6f}\n")
    f.write(f"_cell_angle_alpha {alpha:.6f}\n")
    f.write(f"_cell_angle_beta  {beta:.6f}\n")
    f.write(f"_cell_angle_gamma {gamma:.6f}\n\n")
    f.write("loop_\n")
    f.write("_atom_site_label\n")
    f.write("_atom_site_type_symbol\n")
    f.write("_atom_site_fract_x\n")
    f.write("_atom_site_fract_y\n")
    f.write("_atom_site_fract_z\n")

    for idx, (sym, pos) in enumerate(zip(all_symbols, all_positions), start=1):
        frac = pos @ inv_cell
        f.write(f"{sym}{idx} {sym} {frac[0]:.6f} {frac[1]:.6f} {frac[2]:.6f}\n")

print("\n✅ DONE. Saved to:", OUTPUT_FILE)
