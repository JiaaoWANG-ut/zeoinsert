#!/usr/bin/env python3
"""Probe zeolite/framework pores and export blocked region maps."""

from pore_accessibility import PoreGrid, load_framework_ovito

# ========== 用户参数 ==========
FRAMEWORK_FILE = "frameworks/zeolites/FAU.cif"

# 拓扑封闭孔（完美晶体应 ≈0 blocked）
OUTPUT_BLOCKED_TOPO = "frameworks/zeolites/FAU.blocked_topo.npz"
PROBE_RADIUS_TOPO = 1.2

#  steric 不可达孔（探针 ≈ 客体分子半径，β笼等 confined 区域）
OUTPUT_BLOCKED_STERIC = "frameworks/zeolites/FAU.blocked.npz"
PROBE_RADIUS_STERIC = 2.5

OUTPUT_REPORT = "frameworks/zeolites/FAU.pore_report.txt"

GRID_N = 64
CONNECTIVITY = 26          # 26-邻 flood-fill，减少网格断联 artifact
MIN_CLUSTER_SIZE = 4       # 小于此体素数的 blocked 团簇视为 speckle 并剔除

MANUAL_BLOCKED_BOXES = []


def main():
    print(f"[INFO] Loading framework: {FRAMEWORK_FILE}")
    pos_fw, cell, inv_cell, _ = load_framework_ovito(FRAMEWORK_FILE)
    print(f"[INFO] Framework atoms: {len(pos_fw)}\n")

    configs = [
        ("topological", PROBE_RADIUS_TOPO, OUTPUT_BLOCKED_TOPO),
        ("steric", PROBE_RADIUS_STERIC, OUTPUT_BLOCKED_STERIC),
    ]

    lines = [f"framework: {FRAMEWORK_FILE}", f"grid_n: {GRID_N}",
             f"connectivity: {CONNECTIVITY}", f"min_cluster_size: {MIN_CLUSTER_SIZE}", ""]

    for mode, probe, out_path in configs:
        print(f"===== {mode} (probe={probe} A) =====")
        grid = PoreGrid.build(
            pos_fw, cell, inv_cell,
            grid_n=GRID_N,
            probe_radius=probe,
            manual_boxes=MANUAL_BLOCKED_BOXES,
            framework_file=FRAMEWORK_FILE,
            connectivity=CONNECTIVITY,
            min_cluster_size=MIN_CLUSTER_SIZE,
            mode=mode,
        )
        grid.print_report()
        grid.save(out_path)
        print(f"[INFO] Saved: {out_path}\n")

        lines.append(f"[{mode}] probe={probe} -> {out_path}")
        for k, v in grid.stats.items():
            lines.append(f"  {k}: {v}")
        lines.append("")

    with open(OUTPUT_REPORT, "w") as f:
        f.write("\n".join(lines))
    print(f"[INFO] Saved report: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()
