#!/usr/bin/env python3
import numpy as np
from ovito.io import import_file

from pore_accessibility import PoreGrid

# ========== 用户参数 ==========
FRAMEWORK_FILE = "frameworks/zeolites/FAU.cif"
MOLECULE_FILE  = "molecules/CO2.xyz"
OUTPUT_FILE    = "FAU_mc_all_sa.cif"

N_MOLECULES      = 8      # 同时放置的分子数量

CUTOFF_FW        = 2.2    # 分子 - 框架 最小距离 (Å)
CUTOFF_MM        = 2.2    # 分子 - 分子 最小距离 (Å)

# ========== 封闭孔 / 可达性 ==========
USE_PORE_BLOCKING = True
BLOCKED_FILE = "frameworks/zeolites/FAU.blocked.npz"  # None = 运行时重新探测
GRID_N = 64
PROBE_RADIUS = 1.8
MANUAL_BLOCKED_BOXES = []
BLOCKED_PENALTY = 100.0

MAX_ITERS        = 30000  # 全局 MC 迭代上限
PRINT_EVERY      = 200

# 小步 move 大小
STEP_TRANSL      = 0.8    # Å
STEP_ROT_DEG     = 20.0   # deg

# 模拟退火温度
T0               = 1.0    # 初始温度
T_MIN            = 0.02   # 最小温度

# move 类型概率
P_SMALL_MOVE     = 0.70   # 小平移+小旋转
P_BIG_JUMP       = 0.1   # 全随机新位置+新朝向
P_BIG_ROT        = 0.10   # 仅大角度旋转

# 能量收敛阈值
E_TOL            = 1e-6

# ========== 工具函数 ==========
def random_unit_vector():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)

def random_small_rotation(max_angle_deg):
    theta = np.deg2rad(max_angle_deg) * (2*np.random.rand() - 1.0)
    axis = random_unit_vector()
    kx, ky, kz = axis
    K = np.array([[ 0,  -kz,  ky],
                  [ kz,   0,  -kx],
                  [-ky,  kx,   0]])
    I = np.eye(3)
    R = I + np.sin(theta)*K + (1-np.cos(theta))*(K @ K)
    return R

def random_rotation():
    """ 完全随机刚体旋转（四元数法） """
    u1, u2, u3 = np.random.rand(3)
    q = np.array([
        np.sqrt(1-u1)*np.sin(2*np.pi*u2),
        np.sqrt(1-u1)*np.cos(2*np.pi*u2),
        np.sqrt(u1)*np.sin(2*np.pi*u3),
        np.sqrt(u1)*np.cos(2*np.pi*u3)
    ])
    w, x, y, z = q
    return np.array([
        [1-2*y*y-2*z*z,   2*x*y-2*z*w,   2*x*z+2*y*w],
        [  2*x*y+2*z*w, 1-2*x*x-2*z*z,   2*y*z-2*x*w],
        [  2*x*z-2*y*w,   2*y*z+2*x*w, 1-2*x*x-2*y*y]
    ])

def get_symbols(data):
    if "Element" in data.particles:
        arr = data.particles["Element"].array
        return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]
    elif "Particle Type" in data.particles:
        type_ids = data.particles["Particle Type"].array
        types = data.particles.particle_types.types
        type_map = {t.id: t.name for t in types}
        return [type_map[i] for i in type_ids]
    else:
        raise RuntimeError("No Element info found in data.")

def minimum_image(vecs, cell, inv_cell):
    """ 最小镜像：vecs 可以是 (...,3) """
    frac = vecs @ inv_cell
    frac -= np.round(frac)
    return frac @ cell

def min_distance(A, B, cell, inv_cell):
    """ A:(n,3), B:(m,3) 的最小 PBC 距离 """
    diff = A[:,None,:] - B[None,:,:]
    diff = minimum_image(diff, cell, inv_cell)
    dist = np.linalg.norm(diff, axis=2)
    return float(np.min(dist))

def wrap_to_cell(pos, cell, inv_cell):
    frac = pos @ inv_cell
    frac -= np.floor(frac)
    return frac @ cell

# ========== 读取框架 ==========
print(f"[INFO] Loading framework: {FRAMEWORK_FILE}")
fw = import_file(FRAMEWORK_FILE).compute()
pos_fw = np.array(fw.particles.positions)
sym_fw = get_symbols(fw)

cell = np.array(fw.cell.matrix)[:3,:3]
inv_cell = np.linalg.inv(cell)

print(f"[INFO] Framework atoms: {len(pos_fw)}")

pore_grid = None
if USE_PORE_BLOCKING:
    if BLOCKED_FILE:
        print(f"[INFO] Loading blocked pore map: {BLOCKED_FILE}")
        pore_grid = PoreGrid.load(BLOCKED_FILE)
    else:
        print("[INFO] Probing pore accessibility (no BLOCKED_FILE set)...")
        pore_grid = PoreGrid.build(
            pos_fw, cell, inv_cell,
            grid_n=GRID_N,
            probe_radius=PROBE_RADIUS,
            manual_boxes=MANUAL_BLOCKED_BOXES,
            framework_file=FRAMEWORK_FILE,
        )
    pore_grid.print_report()

def random_center():
    if pore_grid is None:
        return np.random.rand(3)
    return pore_grid.random_accessible_center()

# ========== 读取分子 ==========
print(f"\n[INFO] Loading molecule: {MOLECULE_FILE}")
mol = import_file(MOLECULE_FILE).compute()
pos_mol = np.array(mol.particles.positions)
sym_mol = get_symbols(mol)

# 把分子移到质心
pos_mol -= pos_mol.mean(axis=0)
n_m_atoms = len(pos_mol)

print(f"[INFO] Molecule atoms: {n_m_atoms}")

# ========== 一次性随机初始化所有分子 ==========
centers_frac = np.array([random_center() for _ in range(N_MOLECULES)])
rot_mats     = np.array([random_rotation() for _ in range(N_MOLECULES)])

def world_positions(centers_frac, rot_mats):
    """ 返回 list[N_MOLECULES]，每个是 (n_m_atoms,3) """
    positions = []
    for i in range(N_MOLECULES):
        center_cart = centers_frac[i] @ cell
        R = rot_mats[i]
        positions.append(pos_mol @ R.T + center_cart)
    return positions

def compute_conflicts_and_energy(world_list, centers_frac_local):
    """
    返回 per-molecule 冲突能 conflicts[i] 和总能量 E_tot
    E = sum_over_pairs (max(0, cutoff - d))^2
    """
    conflicts = np.zeros(N_MOLECULES)
    E_tot = 0.0

    if pore_grid is not None:
        for i in range(N_MOLECULES):
            if pore_grid.is_center_blocked(centers_frac_local[i]):
                E_tot += BLOCKED_PENALTY
                conflicts[i] += BLOCKED_PENALTY

    # 分子 - 框架
    for i in range(N_MOLECULES):
        d = min_distance(world_list[i], pos_fw, cell, inv_cell)
        if d < CUTOFF_FW:
            e = (CUTOFF_FW - d)**2
            E_tot += e
            conflicts[i] += e

    # 分子 - 分子（不同分子，不含自身）
    for i in range(N_MOLECULES):
        for j in range(i+1, N_MOLECULES):
            d = min_distance(world_list[i], world_list[j], cell, inv_cell)
            if d < CUTOFF_MM:
                e = (CUTOFF_MM - d)**2
                E_tot += e
                conflicts[i] += e
                conflicts[j] += e

    return conflicts, E_tot

# ========== 模拟退火 MC ==========
print("\n[INFO] Starting global MC with simulated annealing...")

world_list = world_positions(centers_frac, rot_mats)
conflicts, E = compute_conflicts_and_energy(world_list, centers_frac)

print(f"[INFO] Initial E = {E:.4f}")

for it in range(1, MAX_ITERS+1):

    # 退火温度：线性降温 + 地板
    T = max(T_MIN, T0 * (1.0 - it / MAX_ITERS))

    # 若已收敛
    if E < E_TOL:
        print(f"\n✅ Converged at iter {it}, E={E:.6e}")
        break

    # 选冲突最大的分子
    worst = int(np.argmax(conflicts))

    # 记录旧状态
    frac_old = centers_frac[worst].copy()
    R_old    = rot_mats[worst].copy()

    # 决定这一步用哪种 move
    r = np.random.rand()
    if r < P_SMALL_MOVE:
        # 小平移 + 小旋转
        d_cart = STEP_TRANSL * (2*np.random.rand()-1) * random_unit_vector()
        d_frac = d_cart @ inv_cell
        centers_frac[worst] = (frac_old + d_frac) % 1.0
        rot_mats[worst] = random_small_rotation(STEP_ROT_DEG) @ R_old
    elif r < P_SMALL_MOVE + P_BIG_JUMP:
        # 大跳跃：完全随机新位置 + 新旋转
        centers_frac[worst] = random_center()
        rot_mats[worst] = random_rotation()
    else:
        # 仅大角度随机旋转（原地转动）
        centers_frac[worst] = frac_old
        rot_mats[worst] = random_small_rotation(60.0) @ R_old

    # 计算新能量
    world_new = world_positions(centers_frac, rot_mats)
    conflicts_new, E_new = compute_conflicts_and_energy(world_new, centers_frac)

    dE = E_new - E

    # Metropolis 接受准则
    if dE <= 0.0:
        # 接受
        world_list = world_new
        conflicts  = conflicts_new
        E          = E_new
    else:
        p = np.exp(-dE / T)
        if np.random.rand() < p:
            # 以一定概率接受“变差”的 move（跳出 local minima）
            world_list = world_new
            conflicts  = conflicts_new
            E          = E_new
        else:
            # 拒绝，回滚
            centers_frac[worst] = frac_old
            rot_mats[worst]     = R_old

    if it % PRINT_EVERY == 0:
        worst_c = conflicts[worst]
        print(f"[ITER {it}] E={E:.4f}, T={T:.3f}, worst_mol={worst}, conflict={worst_c:.4f}")

else:
    print(f"\n⚠ Reached MAX_ITERS={MAX_ITERS}, E={E:.4f} (may still have mild overlaps).")

# ========== 拼接最终结构 ==========
# world_list 已经是当前构型
all_positions = pos_fw.copy()
all_symbols   = list(sym_fw)

for i in range(N_MOLECULES):
    all_positions = np.vstack((all_positions, world_list[i]))
    all_symbols   += sym_mol

all_positions = wrap_to_cell(all_positions, cell, inv_cell)

# ========== 写 CIF ==========
print("\n[INFO] Writing final CIF:", OUTPUT_FILE)

a = np.linalg.norm(cell[0])
b = np.linalg.norm(cell[1])
c = np.linalg.norm(cell[2])
alpha = np.degrees(np.arccos(np.dot(cell[1],cell[2])/(b*c)))
beta  = np.degrees(np.arccos(np.dot(cell[0],cell[2])/(a*c)))
gamma = np.degrees(np.arccos(np.dot(cell[0],cell[1])/(a*b)))

with open(OUTPUT_FILE, "w") as f:
    f.write("data_mc_all_sa\n")
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
        f.write(
            f"{sym}{idx} {sym} "
            f"{frac[0]:.6f} {frac[1]:.6f} {frac[2]:.6f}\n"
        )

print("\n✅ DONE. Saved to:", OUTPUT_FILE)
