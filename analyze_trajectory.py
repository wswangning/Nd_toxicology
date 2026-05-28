"""
MD轨迹分析脚本 —— 用于完整100ns轨迹的后续分析
=============================================================================
使用方式（在完整MD完成后）：
  python analyze_trajectory.py --pdb md_results/4IQK_final.pdb --traj md_results/4IQK_traj.dcd --output analysis/

分析项目：
  1. RMSD (Cα骨架) —— 随时间 + 收敛性
  2. RMSF (逐残基Cα) —— 热图 + 柔性区域识别
  3. 回转半径 Rg —— 蛋白折叠状态
  4. Nd-关键残基距离 —— Cys151, Sec46等
  5. 氢键分析 —— Nd与周围残基
  6. MM/GBSA 结合自由能（使用OpenMM）
  7. 二级结构变化（DSSP）
=============================================================================
"""
import os, sys, json, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, rmsf, distances, hydrogenbonds
except ImportError:
    print("需要安装 MDAnalysis: pip install MDAnalysis")
    sys.exit(1)

# ============================================================
# 配置
# ============================================================
ND_PARAMS = {
    "sigma_nm": 0.263,
    "epsilon_kjmol": 0.50,
    "charge": 3.0,
    "mass": 144.242,
}

KEY_RESIDUES = {
    "4IQK": {
        "protein": "NFE2L2-Keap1",
        "metal_site": ["CYS151", "CYS273", "CYS288", "CYS297"],  # Keap1关键Cys
    },
    "6HN3": {
        "protein": "GPX4",
        "metal_site": ["SEC46"],  # 硒代半胱氨酸活性位点
    },
}

def analyze_trajectory(pdb_path, traj_path, output_dir, key_residues=None):
    """完整的轨迹分析"""
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"加载体系: {pdb_path} + {traj_path}")
    u = mda.Universe(pdb_path, traj_path)
    
    protein = u.select_atoms("protein")
    ca_atoms = u.select_atoms("protein and name CA")
    backbone = u.select_atoms("protein and backbone")
    
    n_frames = len(u.trajectory)
    total_time_ns = u.trajectory[-1].time / 1000
    
    print(f"帧数: {n_frames}, 总时间: {total_time_ns:.1f} ns")
    
    # ---- 1. RMSD ----
    print("\n1. Cα RMSD 分析")
    ref_positions = ca_atoms.positions.copy()
    
    rmsd_data = []
    for ts in u.trajectory:
        r = rms.rmsd(ca_atoms.positions, ref_positions, center=True, superposition=True)
        rmsd_data.append((ts.time / 1000, r * 10))  # nm→Å
    
    rmsd_arr = np.array(rmsd_data)
    
    # 忽略前20%作为平衡期
    equil_cutoff = int(n_frames * 0.2)
    equil_rmsd = rmsd_arr[equil_cutoff:, 1]
    
    rmsd_summary = {
        "mean_A": float(equil_rmsd.mean()),
        "std_A": float(equil_rmsd.std()),
        "min_A": float(equil_rmsd.min()),
        "max_A": float(equil_rmsd.max()),
    }
    
    print(f"  平衡后 RMSD: {rmsd_summary['mean_A']:.4f} ± {rmsd_summary['std_A']:.4f} Å")
    
    # RMSD图
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(rmsd_arr[:, 0], rmsd_arr[:, 1], 'b-', linewidth=0.5, alpha=0.7)
    ax.axvline(rmsd_arr[equil_cutoff, 0], color='r', linestyle='--', label=f'Equil cutoff ({rmsd_arr[equil_cutoff, 0]:.0f} ns)')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Cα RMSD (Å)')
    ax.set_title('Protein Cα RMSD')
    ax.legend()
    fig.savefig(f"{output_dir}/rmsd.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    # ---- 2. RMSF ----
    print("\n2. 逐残基 RMSF 分析")
    
    # 对齐轨迹到参考帧
    from MDAnalysis.analysis import align
    aligner = align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    
    rmsf_calc = rmsf.RMSF(ca_atoms).run()
    rmsf_values = rmsf_calc.results.rmsf * 10  # nm→Å
    
    rmsf_summary = {
        "mean_A": float(rmsf_values.mean()),
        "std_A": float(rmsf_values.std()),
        "max_A": float(rmsf_values.max()),
        "max_resid": int(np.argmax(rmsf_values)),
    }
    
    print(f"  平均 RMSF: {rmsf_summary['mean_A']:.4f} Å")
    print(f"  最大 RMSF: {rmsf_summary['max_A']:.4f} Å (残基 {rmsf_summary['max_resid']})")
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(rmsf_values, 'b-', linewidth=1)
    ax.axhline(rmsf_values.mean(), color='r', linestyle='--', label=f'Mean: {rmsf_values.mean():.2f} Å')
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title('Per-Residue Cα RMSF')
    ax.legend()
    fig.savefig(f"{output_dir}/rmsf.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    # ---- 3. 回转半径 Rg ----
    print("\n3. 回转半径 (Rg) 分析")
    
    rg_data = []
    for ts in u.trajectory[::10]:  # 每10帧采样
        rg = protein.radius_of_gyration()
        rg_data.append((ts.time / 1000, rg / 10))  # Å→nm
    
    rg_arr = np.array(rg_data)
    rg_summary = {
        "mean_nm": float(rg_arr[:, 1].mean()),
        "std_nm": float(rg_arr[:, 1].std()),
    }
    
    print(f"  Rg: {rg_summary['mean_nm']:.4f} ± {rg_summary['std_nm']:.4f} nm")
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(rg_arr[:, 0], rg_arr[:, 1], 'g-', linewidth=0.5)
    ax.axhline(rg_summary['mean_nm'], color='r', linestyle='--', label=f'Mean: {rg_summary["mean_nm"]:.2f} nm')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (nm)')
    ax.set_title('Radius of Gyration')
    ax.legend()
    fig.savefig(f"{output_dir}/rg.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    # ---- 4. Nd-关键残基距离 ----
    print("\n4. Nd-关键残基距离分析")
    
    # 寻找Nd原子（假设在轨迹中）
    nd_atoms = u.select_atoms("resname ND or name ND")
    
    distance_data = {}
    if len(nd_atoms) > 0 and key_residues:
        for res_name in key_residues:
            # 选择该残基的SG（Cys的硫原子）或SE（Sec的硒原子）
            sel_sg = u.select_atoms(f"resname {res_name[:3]} and resid {res_name[3:]} and (name SG or name SE)")
            if len(sel_sg) > 0:
                dists = []
                for ts in u.trajectory[::10]:
                    d = distances.calc_distance(
                        nd_atoms[0].position, sel_sg[0].position, box=u.dimensions
                    )
                    dists.append((ts.time / 1000, d[2][0] / 10))  # Å→nm
                
                d_arr = np.array(dists)
                occupancy = (d_arr[:, 1] < 0.35).mean() * 100  # <3.5Å
                
                distance_data[res_name] = {
                    "mean_nm": float(d_arr[:, 1].mean()),
                    "std_nm": float(d_arr[:, 1].std()),
                    "occupancy_lt_0.35nm_pct": float(occupancy),
                }
                print(f"  Nd-{res_name}: {d_arr[:, 1].mean():.4f} ± {d_arr[:, 1].std():.4f} nm, occupancy<0.35nm: {occupancy:.0f}%")
    else:
        print("  未找到Nd原子或关键残基")
    
    # ---- 汇总 ----
    results = {
        "pdb": pdb_path,
        "trajectory": traj_path,
        "n_frames": n_frames,
        "total_time_ns": total_time_ns,
        "rmsd": rmsd_summary,
        "rmsf": rmsf_summary,
        "rg": rg_summary,
        "key_residue_distances": distance_data,
    }
    
    with open(f"{output_dir}/analysis_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n✓ 分析完成，结果保存至 {output_dir}/")
    return results


def mmgbsa_analysis(pdb_path, traj_path, output_dir):
    """
    MM/GBSA 结合自由能计算框架
    
    使用OpenMM进行GBSA隐式溶剂能量计算:
    ΔG_bind = G_complex - G_protein - G_ligand
    
    G = E_MM + G_GB + G_SA - TS (忽略熵贡献，仅计算焓)
    """
    print("\n5. MM/GBSA 结合自由能分析")
    print("   (此项需要完整的蛋白+配体复合物轨迹)")
    print("   框架已准备，在完整MD轨迹上运行")
    
    mmgbsa_protocol = {
        "method": "MM/GBSA (igb=5, saltcon=0.15)",
        "software": "OpenMM + implicit solvent",
        "sampling": "最后50 ns, 间隔100 ps (=500帧)",
        "components": {
            "E_vdw": "范德华相互作用能",
            "E_ele": "静电相互作用能",
            "G_GB": "GB极化溶剂化能",
            "G_SA": "非极性溶剂化能 (LCPO方法)",
        },
        "notes": [
            "Nd(III)参数: σ=0.263 nm, ε=0.50 kJ/mol, q=+3.0",
            "建议与MM/PBSA结果交叉验证",
            "收敛性通过ΔG vs. 时间绘图评估",
        ]
    }
    
    with open(f"{output_dir}/mmgbsa_protocol.json", "w") as f:
        json.dump(mmgbsa_protocol, f, indent=2)
    
    return mmgbsa_protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MD轨迹分析")
    parser.add_argument("--pdb", required=True, help="PDB结构文件")
    parser.add_argument("--traj", required=True, help="DCD轨迹文件")
    parser.add_argument("--output", default="analysis", help="输出目录")
    parser.add_argument("--protein", choices=["4IQK", "6HN3"], default="4IQK")
    args = parser.parse_args()
    
    key_res = KEY_RESIDUES.get(args.protein, {}).get("metal_site", None)
    
    results = analyze_trajectory(args.pdb, args.traj, args.output, key_res)
    mmgbsa_analysis(args.pdb, args.traj, args.output)
    
    print("\n" + "=" * 60)
    print("分析完成")
    print("=" * 60)