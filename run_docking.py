"""
Nd(III) 分子对接 —— 基于物理的网格搜索方法
==============================================
针对单原子金属离子的特殊对接方案：
- 在蛋白表面构建3D网格
- 计算Nd(III)在每个格点的相互作用能（LJ + 库仑）
- 选出能量最低的结合姿势
- 生成论文级对接图

Nd(III)参数: σ=0.263 nm, ε=0.50 kJ/mol, q=+3.0
"""
import os, sys, json, time
import numpy as np

os.chdir(r"C:\Users\wangning\AppData\Roaming\Tencent\Marvis\User\oAN1i2eqkIJ6qD9R98ixvZQsjxjI\workspace\conv_19e5eb9a626_881806416319\temp")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 论文级图表样式
rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.linewidth': 1.2,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

from Bio.PDB import PDBParser
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter

# ============================================================
# Nd(III) 参数
# ============================================================
ND_PARAMS = {
    "q": 3.0,           # 电荷 (e)
    "sigma_nm": 0.263,  # LJ σ (nm)
    "epsilon_kjmol": 0.50,  # LJ ε (kJ/mol)
}

# 转换为OpenMM/计算单位
EPSILON_KCAL = ND_PARAMS["epsilon_kjmol"] / 4.184  # kJ/mol → kcal/mol
SIGMA_A = ND_PARAMS["sigma_nm"] * 10  # nm → Å

# 蛋白原子参数（简化的CHARMM/AMBER类型映射）
# 使用平均vdW参数
ATOM_PARAMS = {
    "C":  {"sigma_A": 3.40, "epsilon_kcal": 0.086, "charge": 0.0},   # sp2碳
    "CA": {"sigma_A": 3.40, "epsilon_kcal": 0.086, "charge": 0.0},
    "N":  {"sigma_A": 3.25, "epsilon_kcal": 0.170, "charge": 0.0},   # 酰胺氮
    "O":  {"sigma_A": 2.96, "epsilon_kcal": 0.210, "charge": 0.0},   # 羰基氧
    "S":  {"sigma_A": 3.56, "epsilon_kcal": 0.250, "charge": 0.0},   # 硫(Cys)
    "H":  {"sigma_A": 1.00, "epsilon_kcal": 0.016, "charge": 0.0},
}

def estimate_atom_charge(atom_name, res_name):
    """估算蛋白原子部分电荷"""
    if atom_name == "SG" and res_name == "CYS":
        return -0.23  # 半胱氨酸硫醇
    elif atom_name in ["OD1", "OD2"] and res_name == "ASP":
        return -0.55
    elif atom_name in ["OE1", "OE2"] and res_name == "GLU":
        return -0.55
    elif atom_name in ["NZ"] and res_name == "LYS":
        return 0.31
    elif atom_name in ["NH1", "NH2"] and res_name == "ARG":
        return 0.40
    elif atom_name == "O":
        return -0.55  # 羰基氧
    elif atom_name == "N":
        return -0.20  # 酰胺氮
    else:
        return 0.0

# ============================================================
# Step 1: 读取蛋白结构
# ============================================================
print("=" * 60)
print("Step 1: 读取蛋白结构")
print("=" * 60)

TARGETS = {
    "4IQK": {"name": "Keap1-Nrf2 (NFE2L2)", "pdb": "pdb_structures/4IQK_prepared.pdb"},
    "6HN3": {"name": "GPX4", "pdb": "pdb_structures/6HN3_prepared.pdb"},
}

protein_data = {}

for pdb_id, info in TARGETS.items():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, info["pdb"])
    
    # 提取所有重原子坐标和参数
    coords = []
    atom_info = []  # (sigma, epsilon, charge, res_name, res_id, atom_name)
    
    for atom in structure.get_atoms():
        if atom.element == 'H':
            continue
        name = atom.get_full_id()
        atom_name = name[4][0]
        res_name = name[3][0]
        res_id = name[3][1]
        
        elem = atom.element
        params = ATOM_PARAMS.get(elem, {"sigma_A": 3.00, "epsilon_kcal": 0.10})
        charge = estimate_atom_charge(atom_name, res_name)
        
        coords.append(atom.get_coord())
        atom_info.append({
            "sigma": params["sigma_A"],
            "epsilon": params["epsilon_kcal"],
            "charge": charge,
            "res_name": res_name,
            "res_id": res_id,
            "atom_name": atom_name,
            "element": elem,
        })
    
    coords = np.array(coords)
    bmin = coords.min(axis=0)
    bmax = coords.max(axis=0)
    center = coords.mean(axis=0)
    
    protein_data[pdb_id] = {
        "name": info["name"],
        "coords": coords,
        "atoms": atom_info,
        "n_atoms": len(coords),
        "bbox_min": bmin,
        "bbox_max": bmax,
        "center": center,
    }
    
    print(f"\n{pdb_id} ({info['name']}):")
    print(f"  重原子数: {len(coords)}")
    print(f"  范围: [{bmin[0]:.1f}, {bmin[1]:.1f}, {bmin[2]:.1f}] → [{bmax[0]:.1f}, {bmax[1]:.1f}, {bmax[2]:.1f}] Å")

# ============================================================
# Step 2: 网格对接
# ============================================================
print("\n" + "=" * 60)
print("Step 2: 分子对接 (基于物理的网格搜索)")
print("=" * 60)

def lj_energy(sigma_ij, epsilon_ij, r):
    """Lennard-Jones 12-6 势能"""
    if r < 0.01:
        return 1e10
    sr = sigma_ij / r
    sr6 = sr ** 6
    sr12 = sr6 ** 2
    return epsilon_ij * (sr12 - 2 * sr6)  # 使用2*sr6而非4*(sr12-sr6)因为CHARMM惯例

def coulomb_energy(q1, q2, r, dielectric=4.0):
    """库仑势能 (kcal/mol)，使用距离依赖介电常数"""
    if r < 0.01:
        return 1e10
    # E = 332 * q1 * q2 / (dielectric * r)  (kcal/mol, r in Å)
    return 332.06 * q1 * q2 / (dielectric * r)

def docking_grid_search(protein_coords, atom_info, grid_spacing=1.0, padding=5.0):
    """在蛋白周围的3D网格中搜索Nd(III)最佳结合位置"""
    bmin = protein_coords.min(axis=0) - padding
    bmax = protein_coords.max(axis=0) + padding
    
    # 创建网格
    x = np.arange(bmin[0], bmax[0] + grid_spacing, grid_spacing)
    y = np.arange(bmin[1], bmax[1] + grid_spacing, grid_spacing)
    z = np.arange(bmin[2], bmax[2] + grid_spacing, grid_spacing)
    
    print(f"  网格: {len(x)}×{len(y)}×{len(z)} = {len(x)*len(y)*len(z)} 个格点")
    print(f"  间距: {grid_spacing} Å")
    print(f"  范围: [{bmin[0]:.0f},{bmin[1]:.0f},{bmin[2]:.0f}] → [{bmax[0]:.0f},{bmax[1]:.0f},{bmax[2]:.0f}] Å")
    
    # 构建KD-tree加速
    kdtree = cKDTree(protein_coords)
    
    # Nd参数
    nd_sigma = SIGMA_A
    nd_epsilon = EPSILON_KCAL
    nd_charge = ND_PARAMS["q"]
    
    best_energy = float('inf')
    best_pos = None
    best_details = {}
    
    all_energies = []
    
    # 网格搜索
    n_total = len(x) * len(y) * len(z)
    n_checked = 0
    
    t0 = time.time()
    
    for xi, xv in enumerate(x):
        for yi, yv in enumerate(y):
            for zi, zv in enumerate(z):
                pos = np.array([xv, yv, zv])
                
                # 使用KD-tree找到最近的蛋白原子（判断是否在蛋白内部）
                dist, idx = kdtree.query(pos, k=1)
                
                # 跳过蛋白内部（距离任何原子<2.0Å）
                if dist < 2.0:
                    continue
                
                # 找到所有在截断距离内的原子
                nearby_idx = kdtree.query_ball_point(pos, r=12.0)
                
                if len(nearby_idx) == 0:
                    continue
                
                e_total = 0.0
                e_vdw = 0.0
                e_ele = 0.0
                
                for idx_a in nearby_idx:
                    atom = atom_info[idx_a]
                    r = np.linalg.norm(pos - protein_coords[idx_a])
                    
                    if r < 0.5:
                        e_total = 1e10
                        break
                    
                    # LJ能量
                    sigma_ij = (nd_sigma + atom["sigma"]) / 2.0
                    epsilon_ij = np.sqrt(nd_epsilon * atom["epsilon"])
                    e_vdw += lj_energy(sigma_ij, epsilon_ij, r)
                    
                    # 库仑能量
                    if abs(atom["charge"]) > 0.001:
                        e_ele += coulomb_energy(nd_charge, atom["charge"], r)
                
                e_total = e_vdw + e_ele
                all_energies.append((e_total, pos.copy(), e_vdw, e_ele))
                
                if e_total < best_energy:
                    best_energy = e_total
                    best_pos = pos.copy()
                    best_details = {
                        "E_total": e_total,
                        "E_vdw": e_vdw,
                        "E_ele": e_ele,
                    }
                
                n_checked += 1
                if n_checked % 500000 == 0:
                    elapsed = time.time() - t0
                    progress = n_checked / n_total * 100
                    print(f"  进度: {progress:.1f}% ({n_checked}/{n_total}), 耗时: {elapsed:.1f}s, 最优E: {best_energy:.2f}")
    
    elapsed = time.time() - t0
    print(f"\n  完成! 检查了 {n_checked} 个格点, 耗时 {elapsed:.1f}s")
    
    # 排序并取top姿势
    all_energies.sort(key=lambda x: x[0])
    return all_energies[:50], best_pos, best_details

# 对每个蛋白执行对接
docking_results = {}

for pdb_id, data in protein_data.items():
    print(f"\n--- {pdb_id} ({data['name']}) ---")
    
    # 先用粗网格快速扫描
    print("  粗网格扫描 (2.0 Å)...")
    coarse_results, _, _ = docking_grid_search(
        data["coords"], data["atoms"], grid_spacing=2.0, padding=8.0
    )
    
    # 在最优位置附近用细网格精细搜索
    best_coarse_pos = coarse_results[0][1]
    print(f"\n  粗网格最优位置: {best_coarse_pos}")
    print(f"  粗网格最优能量: E_total={coarse_results[0][0]:.2f}, E_vdw={coarse_results[0][2]:.2f}, E_ele={coarse_results[0][3]:.2f} kcal/mol")
    
    # 精细搜索（在最优位置±10Å范围内，0.5Å间距）
    local_center = best_coarse_pos
    half_range = 10.0
    
    mask = np.all(np.abs(data["coords"] - local_center) < half_range * 2, axis=1)
    local_coords = data["coords"][mask]
    
    # 重建原子索引映射
    local_indices = np.where(mask)[0]
    local_atoms = [data["atoms"][i] for i in local_indices]
    
    x = np.arange(local_center[0] - half_range, local_center[0] + half_range, 0.5)
    y = np.arange(local_center[1] - half_range, local_center[1] + half_range, 0.5)
    z = np.arange(local_center[2] - half_range, local_center[2] + half_range, 0.5)
    
    print(f"\n  精细搜索 (0.5 Å, {len(x)}×{len(y)}×{len(z)} 格点)...")
    
    nd_sigma = SIGMA_A
    nd_epsilon = EPSILON_KCAL
    nd_charge = ND_PARAMS["q"]
    
    fine_results = []
    best_fine_e = float('inf')
    best_fine_pos = None
    
    for xv in x:
        for yv in y:
            for zv in z:
                pos = np.array([xv, yv, zv])
                e_vdw = 0.0
                e_ele = 0.0
                
                for i, coord in enumerate(local_coords):
                    r = np.linalg.norm(pos - coord)
                    if r < 0.5:
                        e_vdw = 1e10
                        break
                    
                    atom = local_atoms[i]
                    sigma_ij = (nd_sigma + atom["sigma"]) / 2.0
                    epsilon_ij = np.sqrt(nd_epsilon * atom["epsilon"])
                    e_vdw += lj_energy(sigma_ij, epsilon_ij, r)
                    
                    if abs(atom["charge"]) > 0.001:
                        e_ele += coulomb_energy(nd_charge, atom["charge"], r)
                
                e_total = e_vdw + e_ele
                
                if e_total < best_fine_e:
                    best_fine_e = e_total
                    best_fine_pos = pos.copy()
    
    # 查找结合位点附近的残基
    kdtree = cKDTree(data["coords"])
    nearby = kdtree.query_ball_point(best_fine_pos, r=5.0)
    
    nearby_residues = set()
    for idx in nearby:
        atom = data["atoms"][idx]
        nearby_residues.add(f"{atom['res_name']}{atom['res_id']}")
    
    # 查找金属结合残基（Cys/Sec/His/Glu/Asp）
    metal_binding = []
    for idx in nearby:
        atom = data["atoms"][idx]
        if atom["res_name"] in ["CYS", "SEC", "HIS", "GLU", "ASP"]:
            r = np.linalg.norm(best_fine_pos - data["coords"][idx])
            if r < 4.0:
                metal_binding.append({
                    "residue": f"{atom['res_name']}{atom['res_id']}",
                    "atom": atom["atom_name"],
                    "distance_A": round(float(r), 2),
                })
    
    result = {
        "pdb_id": pdb_id,
        "name": data["name"],
        "docking_energy_kcal_mol": round(float(best_fine_e), 2),
        "best_position": best_fine_pos.tolist(),
        "coarse_energy": round(float(coarse_results[0][0]), 2),
        "nearby_residues": sorted(nearby_residues),
        "metal_binding_residues": metal_binding,
        "top_10_coarse": [
            {"energy": round(float(e), 2), "pos": p.tolist(), "E_vdw": round(float(v), 2), "E_ele": round(float(el), 2)}
            for e, p, v, el in coarse_results[:10]
        ]
    }
    
    docking_results[pdb_id] = result
    
    print(f"\n  ✓ 最优结合位点: {best_fine_pos}")
    print(f"  ✓ 对接能量: E_total={best_fine_e:.2f} kcal/mol")
    print(f"  ✓ 5Å内残基: {', '.join(sorted(nearby_residues))}")
    if metal_binding:
        print(f"  ✓ 金属结合残基 (<4Å):")
        for mb in metal_binding:
            print(f"      {mb['residue']}/{mb['atom']}: {mb['distance_A']} Å")

# 保存对接结果
os.makedirs("docking_results", exist_ok=True)
with open("docking_results/grid_docking_results.json", "w") as f:
    json.dump(docking_results, f, indent=2, default=str)

print("\n✓ 对接结果已保存: docking_results/grid_docking_results.json")

# ============================================================
# Step 3: 生成对接结果图
# ============================================================
print("\n" + "=" * 60)
print("Step 3: 生成对接结果图（论文发表级别）")
print("=" * 60)

os.makedirs("figures", exist_ok=True)

for pdb_id, result in docking_results.items():
    data = protein_data[pdb_id]
    
    # --- 图1: 结合能分布热图 (2D投影) ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    titles = ['XY Projection', 'XZ Projection', 'YZ Projection']
    projections = [(0,1), (0,2), (1,2)]
    
    for ax, (pi, pj), title in zip(axes, projections, titles):
        # 绘制蛋白原子（简化：骨架CA）
        ca_mask = np.array([a["atom_name"] == "CA" for a in data["atoms"]])
        ca_coords = data["coords"][ca_mask]
        
        ax.scatter(ca_coords[:, pi], ca_coords[:, pj], s=0.3, c='lightgray', alpha=0.5, rasterized=True)
        
        # 绘制最佳对接位置
        best = np.array(result["best_position"])
        ax.scatter(best[pi], best[pj], s=150, c='red', marker='*', 
                  edgecolors='darkred', linewidths=1.5, zorder=5, label='Nd(III) best pose')
        
        # 标注金属结合残基
        for mb in result["metal_binding_residues"]:
            for idx, atom in enumerate(data["atoms"]):
                if atom["res_name"] + str(atom["res_id"]) == mb["residue"] and atom["atom_name"] == "CA":
                    coord = data["coords"][idx]
                    ax.annotate(mb["residue"], (coord[pi], coord[pj]),
                              fontsize=6, color='blue',
                              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                    break
        
        ax.set_xlabel(f'{["X","X","Y"][["XY","XZ","YZ"].index(title.split()[0])]} (Å)')
        ax.set_ylabel(f'{["Y","Z","Z"][["XY","XZ","YZ"].index(title.split()[0])]} (Å)')
        ax.set_title(title, fontweight='bold')
        ax.legend(loc='upper right', fontsize=7)
        ax.set_aspect('equal')
    
    fig.suptitle(f'Nd(III) Docking Pose — {result["name"]} ({pdb_id})\n'
                f'Binding Energy: {result["docking_energy_kcal_mol"]:.2f} kcal/mol',
                fontweight='bold')
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_docking_pose.png", dpi=300)
    plt.close(fig)
    print(f"  ✓ {pdb_id}_docking_pose.png")
    
    # --- 图2: 粗网格能量分布 ---
    coarse_energies = [r["energy"] for r in result["top_10_coarse"]]
    coarse_labels = [f"Pose {i+1}" for i in range(len(coarse_energies))]
    
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = ['#d62728' if i == 0 else '#1f77b4' for i in range(len(coarse_energies))]
    bars = ax.barh(range(len(coarse_energies)), coarse_energies, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(coarse_energies)))
    ax.set_yticklabels(coarse_labels)
    ax.set_xlabel('Binding Energy (kcal/mol)')
    ax.set_title(f'Top 10 Docking Poses — {result["name"]} ({pdb_id})', fontweight='bold')
    ax.axvline(x=coarse_energies[0], color='red', linestyle='--', alpha=0.5)
    ax.invert_yaxis()
    
    for i, (e, bar) in enumerate(zip(coarse_energies, bars)):
        ax.text(e + 0.5, bar.get_y() + bar.get_height()/2, f'{e:.1f}', va='center', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_docking_energies.png", dpi=300)
    plt.close(fig)
    print(f"  ✓ {pdb_id}_docking_energies.png")
    
    # --- 图3: 结合位点残基相互作用图 ---
    if result["metal_binding_residues"]:
        fig, ax = plt.subplots(figsize=(8, 4))
        
        residues = [mb["residue"] for mb in result["metal_binding_residues"]]
        distances = [mb["distance_A"] for mb in result["metal_binding_residues"]]
        
        colors_bar = ['#2ca02c' if d < 3.0 else '#ff7f0e' if d < 3.5 else '#d62728' for d in distances]
        
        bars = ax.barh(range(len(residues)), distances, color=colors_bar, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(residues)))
        ax.set_yticklabels(residues)
        ax.set_xlabel('Distance from Nd(III) (Å)')
        ax.set_title(f'Metal-Binding Site Interactions — {result["name"]} ({pdb_id})', fontweight='bold')
        ax.axvline(x=3.0, color='green', linestyle='--', alpha=0.5, label='Coordination (3.0 Å)')
        ax.axvline(x=3.5, color='red', linestyle='--', alpha=0.5, label='Weak interaction (3.5 Å)')
        ax.legend(fontsize=8)
        ax.invert_yaxis()
        
        for i, (d, bar) in enumerate(zip(distances, bars)):
            ax.text(d + 0.05, bar.get_y() + bar.get_height()/2, f'{d:.2f} Å', va='center', fontsize=9)
        
        plt.tight_layout()
        fig.savefig(f"figures/{pdb_id}_binding_site.png", dpi=300)
        plt.close(fig)
        print(f"  ✓ {pdb_id}_binding_site.png")

# ============================================================
# 最终汇总
# ============================================================
print("\n" + "=" * 60)
print("分子对接完成")
print("=" * 60)

for pdb_id, result in docking_results.items():
    print(f"\n{pdb_id} ({result['name']}):")
    print(f"  对接能量: {result['docking_energy_kcal_mol']:.2f} kcal/mol")
    print(f"  最优位置: {result['best_position']}")
    print(f"  金属结合残基: {[mb['residue'] for mb in result['metal_binding_residues']]}")

print("\n生成的文件:")
print("  docking_results/grid_docking_results.json — 完整对接结果")
for pdb_id in TARGETS:
    print(f"  figures/{pdb_id}_docking_pose.png — 对接姿势图")
    print(f"  figures/{pdb_id}_docking_energies.png — Top10能量图")
    print(f"  figures/{pdb_id}_binding_site.png — 结合位点图")

print("\nOK")