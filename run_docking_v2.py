"""
Nd(III) 分子对接 —— 聚焦金属结合位点的智能搜索
==================================================
策略：
1. 识别蛋白中所有Cys/Sec/His/Glu/Asp（金属结合热点残基）
2. 聚类这些残基找到潜在结合区域
3. 仅在这些热点区域做精细网格搜索
4. 大幅减少计算量（从百万级格点降到千级）
"""
import os, sys, json, time
import numpy as np

os.chdir(r"C:\Users\wangning\AppData\Roaming\Tencent\Marvis\User\oAN1i2eqkIJ6qD9R98ixvZQsjxjI\workspace\conv_19e5eb9a626_881806416319\temp")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({
    'font.family': 'sans-serif', 'font.size': 10,
    'axes.linewidth': 1.2, 'axes.labelsize': 12,
    'axes.titlesize': 13, 'xtick.labelsize': 9,
    'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

from Bio.PDB import PDBParser
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN

# ============================================================
# 参数
# ============================================================
SIGMA_A = 2.63   # Nd(III) sigma in Å (0.263 nm)
EPSILON_KCAL = 0.50 / 4.184  # Nd(III) epsilon in kcal/mol
ND_CHARGE = 3.0

ATOM_PARAMS = {
    "C": {"sigma": 3.40, "epsilon": 0.086},
    "CA": {"sigma": 3.40, "epsilon": 0.086},
    "N": {"sigma": 3.25, "epsilon": 0.170},
    "O": {"sigma": 2.96, "epsilon": 0.210},
    "S": {"sigma": 3.56, "epsilon": 0.250},
    "H": {"sigma": 1.00, "epsilon": 0.016},
    "SE": {"sigma": 3.80, "epsilon": 0.300},
}

def estimate_charge(atom_name, res_name):
    """估算蛋白原子部分电荷"""
    if atom_name == "SG" and res_name == "CYS": return -0.35
    if atom_name == "SEG" and res_name == "SEC": return -0.40
    if atom_name in ["OD1","OD2"] and res_name == "ASP": return -0.55
    if atom_name in ["OE1","OE2"] and res_name == "GLU": return -0.55
    if atom_name in ["ND1","NE2"] and res_name == "HIS": return -0.40
    if atom_name in ["NZ"] and res_name == "LYS": return 0.31
    if atom_name in ["NH1","NH2"] and res_name == "ARG": return 0.40
    if atom_name == "O": return -0.55
    if atom_name == "N": return -0.20
    return 0.0

def lj_energy(sigma_ij, epsilon_ij, r):
    if r < 0.5: return 1e10
    sr = sigma_ij / r
    sr6 = sr ** 6
    return epsilon_ij * (sr6**2 - 2 * sr6)

def coulomb_energy(q1, q2, r, dielectric=4.0):
    if r < 0.5: return 1e10
    return 332.06 * q1 * q2 / (dielectric * r)

# ============================================================
# Step 1: 读取蛋白并识别金属结合热点
# ============================================================
print("=" * 60)
print("Step 1: 读取蛋白结构 & 识别金属结合热点区域")
print("=" * 60)

TARGETS = {
    "4IQK": {"name": "Keap1-Nrf2 (NFE2L2)", "pdb": "pdb_structures/4IQK_prepared.pdb"},
    "6HN3": {"name": "GPX4", "pdb": "pdb_structures/6HN3_prepared.pdb"},
}

METAL_BINDING_RES = {"CYS", "SEC", "HIS", "GLU", "ASP", "MET"}

protein_data = {}

for pdb_id, info in TARGETS.items():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, info["pdb"])
    
    all_coords = []
    all_atoms_info = []
    hot_coords = []  # 金属结合残基重原子坐标
    hot_residues = []  # 热点残基信息
    
    for atom in structure.get_atoms():
        if atom.element == 'H': continue
        atom_name = atom.get_name()
        residue = atom.get_parent()
        res_name = residue.get_resname()
        res_id = residue.get_id()[1]
        elem = atom.element
        
        params = ATOM_PARAMS.get(elem, {"sigma": 3.00, "epsilon": 0.10})
        charge = estimate_charge(atom_name, res_name)
        coord = atom.get_coord()
        
        all_coords.append(coord)
        all_atoms_info.append({
            "sigma": params["sigma"], "epsilon": params["epsilon"],
            "charge": charge, "res_name": res_name,
            "res_id": res_id, "atom_name": atom_name, "element": elem,
        })
        
        if res_name in METAL_BINDING_RES:
            hot_coords.append(coord)
            hot_residues.append({
                "res_name": res_name, "res_id": res_id,
                "atom_name": atom_name, "coord": coord,
            })
    
    all_coords = np.array(all_coords)
    hot_coords = np.array(hot_coords)
    
    # 聚类热点残基找潜在结合位点
    hot_res_set = sorted(set(f"{r['res_name']}{r['res_id']}" for r in hot_residues))
    
    # 每个热点残基的代表坐标（取重原子的中心）
    hot_centers = {}
    for r in hot_residues:
        key = f"{r['res_name']}{r['res_id']}"
        if key not in hot_centers:
            hot_centers[key] = []
        hot_centers[key].append(r["coord"])
    
    hot_keys = list(hot_centers.keys())
    hot_center_coords = np.array([np.mean(hot_centers[k], axis=0) for k in hot_keys])
    
    # DBSCAN聚类（8Å半径内至少2个热点残基为一个cluster）
    if len(hot_center_coords) >= 2:
        clustering = DBSCAN(eps=8.0, min_samples=2).fit(hot_center_coords)
        labels = clustering.labels_
    else:
        labels = np.zeros(len(hot_center_coords), dtype=int)
    
    clusters = {}
    for i, label in enumerate(labels):
        if label == -1: continue
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(hot_keys[i])
    
    protein_data[pdb_id] = {
        "name": info["name"],
        "all_coords": all_coords,
        "all_atoms": all_atoms_info,
        "n_atoms": len(all_coords),
        "hot_residues": hot_res_set,
        "n_hot": len(hot_res_set),
        "clusters": clusters,
        "hot_centers": hot_centers,
        "hot_center_coords": hot_center_coords,
        "hot_keys": hot_keys,
        "labels": labels,
    }
    
    print(f"\n{pdb_id} ({info['name']}):")
    print(f"  总重原子: {len(all_coords)}")
    print(f"  金属结合残基(Cys/Sec/His/Glu/Asp/Met): {len(hot_res_set)}")
    print(f"  识别到 {len(clusters)} 个潜在结合位点聚类:")
    for cid, cres in clusters.items():
        print(f"    Cluster {cid}: {cres}")

# ============================================================
# Step 2: 聚焦热点区域做精细网格搜索
# ============================================================
print("\n" + "=" * 60)
print("Step 2: 聚焦热点区域精细网格搜索 (0.5 Å)")
print("=" * 60)

docking_results = {}

for pdb_id, data in protein_data.items():
    print(f"\n--- {pdb_id} ({data['name']}) ---")
    
    all_coords = data["all_coords"]
    all_atoms = data["all_atoms"]
    
    # 对每个cluster做精细搜索
    all_poses = []
    
    for cid, cluster_res in data["clusters"].items():
        # 计算cluster的几何中心
        cluster_coords = []
        for rkey in cluster_res:
            cluster_coords.extend(data["hot_centers"][rkey])
        cluster_coords = np.array(cluster_coords)
        cluster_center = cluster_coords.mean(axis=0)
        
        print(f"  Cluster {cid}: {cluster_res}, center={cluster_center.round(1)}")
        
        # 在cluster周围做精细搜索 (范围8Å, 间距0.5Å)
        search_range = 8.0
        spacing = 0.5
        
        x = np.arange(cluster_center[0]-search_range, cluster_center[0]+search_range, spacing)
        y = np.arange(cluster_center[1]-search_range, cluster_center[1]+search_range, spacing)
        z = np.arange(cluster_center[2]-search_range, cluster_center[2]+search_range, spacing)
        
        print(f"    搜索网格: {len(x)}×{len(y)}×{len(z)} = {len(x)*len(y)*len(z)} 格点")
        
        # 用KD-tree加速
        kdtree = cKDTree(all_coords)
        
        cluster_best_e = float('inf')
        cluster_best_pos = None
        cluster_best_details = {}
        
        t0 = time.time()
        n_checked = 0
        
        for xv in x:
            for yv in y:
                for zv in z:
                    pos = np.array([xv, yv, zv])
                    
                    # 跳过蛋白内部
                    dist, _ = kdtree.query(pos, k=1)
                    if dist < 2.0: continue
                    
                    # 找12Å内原子
                    nearby = kdtree.query_ball_point(pos, r=12.0)
                    if len(nearby) == 0: continue
                    
                    e_vdw = 0.0
                    e_ele = 0.0
                    
                    for idx in nearby:
                        atom = all_atoms[idx]
                        r = np.linalg.norm(pos - all_coords[idx])
                        if r < 0.5:
                            e_vdw = 1e10
                            break
                        
                        sigma_ij = (SIGMA_A + atom["sigma"]) / 2.0
                        epsilon_ij = np.sqrt(EPSILON_KCAL * atom["epsilon"])
                        e_vdw += lj_energy(sigma_ij, epsilon_ij, r)
                        
                        if abs(atom["charge"]) > 0.001:
                            e_ele += coulomb_energy(ND_CHARGE, atom["charge"], r)
                    
                    e_total = e_vdw + e_ele
                    n_checked += 1
                    
                    if e_total < cluster_best_e:
                        cluster_best_e = e_total
                        cluster_best_pos = pos.copy()
                        cluster_best_details = {
                            "E_total": e_total, "E_vdw": e_vdw, "E_ele": e_ele,
                            "cluster": cid, "cluster_residues": cluster_res,
                        }
                
                if n_checked % 10000 == 0:
                    elapsed = time.time() - t0
                    # Don't print too much
        
        elapsed = time.time() - t0
        print(f"    检查了 {n_checked} 格点, 耗时 {elapsed:.1f}s")
        print(f"    最优能量: E_total={cluster_best_e:.2f}, E_vdw={cluster_best_details['E_vdw']:.2f}, E_ele={cluster_best_details['E_ele']:.2f} kcal/mol")
        
        all_poses.append({
            "energy": cluster_best_e,
            "position": cluster_best_pos.tolist(),
            "details": cluster_best_details,
        })
    
    # 选择全局最优姿势
    all_poses.sort(key=lambda x: x["energy"])
    best_pose = all_poses[0]
    
    # 识别最佳姿势附近的残基
    best_pos = np.array(best_pose["position"])
    kdtree = cKDTree(all_coords)
    nearby_idx = kdtree.query_ball_point(best_pos, r=5.0)
    
    nearby_res = set()
    metal_binding = []
    for idx in nearby_idx:
        atom = all_atoms[idx]
        nearby_res.add(f"{atom['res_name']}{atom['res_id']}")
        if atom["res_name"] in METAL_BINDING_RES:
            r = np.linalg.norm(best_pos - all_coords[idx])
            if r < 4.5:
                metal_binding.append({
                    "residue": f"{atom['res_name']}{atom['res_id']}",
                    "atom": atom["atom_name"],
                    "distance_A": round(float(r), 2),
                })
    
    # 去重金属结合残基（取最近距离）
    mb_dict = {}
    for mb in metal_binding:
        key = mb["residue"]
        if key not in mb_dict or mb["distance_A"] < mb_dict[key]["distance_A"]:
            mb_dict[key] = mb
    metal_binding = sorted(mb_dict.values(), key=lambda x: x["distance_A"])
    
    result = {
        "pdb_id": pdb_id,
        "name": data["name"],
        "docking_energy_kcal_mol": round(float(best_pose["energy"]), 2),
        "best_position": best_pose["position"],
        "best_cluster": best_pose["details"]["cluster"],
        "best_cluster_residues": best_pose["details"]["cluster_residues"],
        "nearby_residues_5A": sorted(nearby_res),
        "metal_binding_residues": metal_binding,
        "all_cluster_poses": all_poses,
        "n_clusters": len(data["clusters"]),
        "n_hot_residues": data["n_hot"],
    }
    
    docking_results[pdb_id] = result
    
    print(f"\n  ✓ 全局最优: {best_pose['position']}")
    print(f"  ✓ 对接能量: E_total={best_pose['energy']:.2f} kcal/mol")
    print(f"  ✓ 最佳聚类: Cluster {best_pose['details']['cluster']}: {best_pose['details']['cluster_residues']}")
    print(f"  ✓ 5Å内残基: {', '.join(sorted(nearby_res)[:20])}")
    if metal_binding:
        print(f"  ✓ 金属结合残基 (<4.5Å):")
        for mb in metal_binding[:8]:
            print(f"      {mb['residue']}/{mb['atom']}: {mb['distance_A']} Å")

os.makedirs("docking_results", exist_ok=True)
with open("docking_results/grid_docking_results.json", "w") as f:
    json.dump(docking_results, f, indent=2, default=str)

print("\n✓ 对接结果已保存")

# ============================================================
# Step 3: 论文级图表
# ============================================================
print("\n" + "=" * 60)
print("Step 3: 生成论文级图表")
print("=" * 60)

os.makedirs("figures", exist_ok=True)

for pdb_id, result in docking_results.items():
    data = protein_data[pdb_id]
    
    # --- 图1: 对接姿势3视图 ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    proj_names = ['XY Projection', 'XZ Projection', 'YZ Projection']
    proj_pairs = [(0,1), (0,2), (1,2)]
    
    best = np.array(result["best_position"])
    
    for ax, (pi, pj), title in zip(axes, proj_pairs, proj_names):
        # 蛋白CA骨架
        ca_mask = np.array([a["atom_name"] == "CA" for a in data["all_atoms"]])
        ca_coords = data["all_coords"][ca_mask]
        ax.scatter(ca_coords[:,pi], ca_coords[:,pj], s=0.2, c='lightgray', alpha=0.4, rasterized=True)
        
        # 热点残基高亮
        hot_mask = np.array([
            a["res_name"] in METAL_BINDING_RES and a["atom_name"] == "CA"
            for a in data["all_atoms"]
        ])
        hot_ca = data["all_coords"][hot_mask]
        ax.scatter(hot_ca[:,pi], hot_ca[:,pj], s=3, c='steelblue', alpha=0.6, label='Metal-binding residues', rasterized=True)
        
        # Nd(III)最佳位置
        ax.scatter(best[pi], best[pj], s=200, c='red', marker='*',
                  edgecolors='darkred', linewidths=2, zorder=10, label='Nd(III)')
        
        # 最佳聚类残基标注
        for rkey in result["best_cluster_residues"][:6]:
            for idx, atom in enumerate(data["all_atoms"]):
                if f"{atom['res_name']}{atom['res_id']}" == rkey and atom["atom_name"] == "CA":
                    c = data["all_coords"][idx]
                    ax.annotate(rkey, (c[pi], c[pj]), fontsize=6, color='darkblue',
                              bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.8))
                    break
        
        ax.set_xlabel(f'{["X","X","Y"][proj_pairs.index((pi,pj))]} (Å)')
        ax.set_ylabel(f'{["Y","Z","Z"][proj_pairs.index((pi,pj))]} (Å)')
        ax.set_title(title, fontweight='bold')
        ax.legend(loc='best', fontsize=7)
        ax.set_aspect('equal')
    
    fig.suptitle(f'Molecular Docking of Nd(III) to {result["name"]} ({pdb_id})\n'
                f'$\\Delta G_{{bind}}$ = {result["docking_energy_kcal_mol"]:.1f} kcal/mol',
                fontweight='bold', fontsize=13)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_docking_pose.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {pdb_id}_docking_pose.png")
    
    # --- 图2: 结合位点详细相互作用 ---
    if result["metal_binding_residues"]:
        fig, ax = plt.subplots(figsize=(9, 5))
        
        residues = [mb["residue"] for mb in result["metal_binding_residues"][:12]]
        distances = [mb["distance_A"] for mb in result["metal_binding_residues"][:12]]
        
        cmap_colors = plt.cm.RdYlGn_r([min(1.0, (d-2.0)/3.0) for d in distances])
        
        bars = ax.barh(range(len(residues)), distances, color=cmap_colors, edgecolor='black', linewidth=0.8)
        ax.set_yticks(range(len(residues)))
        ax.set_yticklabels(residues, fontfamily='monospace')
        ax.set_xlabel('Distance to Nd(III) (Å)', fontsize=12)
        ax.set_title(f'Metal-Binding Site Interaction Distances — {pdb_id}', fontweight='bold', fontsize=13)
        
        ax.axvspan(2.0, 3.0, alpha=0.08, color='green', label='Coordination bond (2.0-3.0 Å)')
        ax.axvspan(3.0, 4.0, alpha=0.08, color='orange', label='Electrostatic (3.0-4.0 Å)')
        ax.axvspan(4.0, 5.0, alpha=0.08, color='red', label='Weak interaction (>4.0 Å)')
        
        ax.legend(loc='lower right', fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, max(distances) + 1)
        
        for i, (d, bar) in enumerate(zip(distances, bars)):
            ax.text(d + 0.08, bar.get_y() + bar.get_height()/2, f'{d:.2f}',
                   va='center', fontsize=9, fontweight='bold')
        
        plt.tight_layout()
        fig.savefig(f"figures/{pdb_id}_binding_interactions.png", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ {pdb_id}_binding_interactions.png")
    
    # --- 图3: 各聚类位点能量对比 ---
    cluster_labels = [f"Cluster {p['details']['cluster']}" for p in result["all_cluster_poses"]]
    cluster_energies = [p["energy"] for p in result["all_cluster_poses"]]
    
    fig, ax = plt.subplots(figsize=(8, 4))
    colors_cl = ['#d62728' if i == 0 else '#1f77b4' for i in range(len(cluster_energies))]
    bars = ax.bar(range(len(cluster_energies)), cluster_energies, color=colors_cl, edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(cluster_energies)))
    ax.set_xticklabels(cluster_labels, rotation=30, fontsize=9)
    ax.set_ylabel('Binding Energy (kcal/mol)', fontsize=12)
    ax.set_title(f'Cluster Site Comparison — {pdb_id}', fontweight='bold', fontsize=13)
    
    for i, (e, bar) in enumerate(zip(cluster_energies, bars)):
        ax.text(i, e + 1, f'{e:.1f}', ha='center', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_cluster_comparison.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {pdb_id}_cluster_comparison.png")

# ============================================================
# 汇总
# ============================================================
print("\n" + "=" * 60)
print("分子对接完成")
print("=" * 60)
for pdb_id, result in docking_results.items():
    print(f"\n{pdb_id}: E_bind={result['docking_energy_kcal_mol']:.2f} kcal/mol")
    print(f"  Cluster: {result['best_cluster_residues']}")
    print(f"  Key residues: {[mb['residue'] for mb in result['metal_binding_residues'][:6]]}")

print("\nOK")
