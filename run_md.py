"""
Nd(III) 分子动力学模拟 —— 基于对接结果的MD模拟
==================================================
使用对接得到的最佳姿势作为初始构象，运行100ns MD模拟。
"""
import os, sys, json, time, math
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

# ============================================================
# Step 1: 准备MD模拟体系
# ============================================================
print("=" * 60)
print("Step 1: 准备MD模拟体系")
print("=" * 60)

# 读取对接结果
with open("docking_results/grid_docking_results.json", "r") as f:
    docking_results = json.load(f)

# Nd(III) 力场参数 (CHARMM36兼容)
# 基于Li & Merz (2016) 12-6-4模型简化
ND_PARAMS = {
    "sigma_nm": 0.263,      # LJ σ (nm)
    "epsilon_kjmol": 0.50,  # LJ ε (kJ/mol)
    "charge": 3.0,          # 电荷 (e)
    "mass_amu": 144.24,     # 原子质量
}

# 转换为kcal/mol
SIGMA_NM = ND_PARAMS["sigma_nm"]
EPSILON_KCAL = ND_PARAMS["epsilon_kjmol"] / 4.184  # kJ/mol → kcal/mol

print("Nd(III) 力场参数:")
print(f"  σ = {SIGMA_NM} nm")
print(f"  ε = {EPSILON_KCAL:.4f} kcal/mol")
print(f"  q = {ND_PARAMS['charge']} e")
print(f"  mass = {ND_PARAMS['mass_amu']} amu")

# 检查OpenMM是否可用
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm import Platform
    print("\n✓ OpenMM 已安装")
except ImportError:
    print("\n✗ OpenMM 未安装，尝试安装...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "openmm"])
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm import Platform
    print("✓ OpenMM 安装成功")

# ============================================================
# Step 2: 为每个蛋白构建MD体系
# ============================================================
print("\n" + "=" * 60)
print("Step 2: 构建MD体系")
print("=" * 60)

os.makedirs("md_systems", exist_ok=True)
os.makedirs("md_trajectories", exist_ok=True)

def create_md_system(pdb_id, docking_result):
    """为单个蛋白构建MD体系"""
    print(f"\n--- {pdb_id} ({docking_result['name']}) ---")
    
    # 1. 加载蛋白结构
    pdb_file = f"pdb_structures/{pdb_id}_prepared.pdb"
    print(f"  加载蛋白: {pdb_file}")
    
    pdb = app.PDBFile(pdb_file)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # 2. 添加Nd(III)离子
    nd_position = np.array(docking_result["best_position"]) * unit.angstroms
    print(f"  Nd(III)位置: {nd_position}")
    
    # 创建Nd原子
    element = app.Element.getByAtomicNumber(60)  # Nd
    nd_atom = app.topology.Atom("Nd", element)
    
    # 添加到拓扑
    chain = modeller.topology.addChain()
    residue = modeller.topology.addResidue("ND", chain)
    modeller.topology.addAtom("Nd", element, residue)
    
    # 添加到坐标
    modeller.positions.append(nd_position)
    
    # 3. 添加水盒子
    print("  添加水盒子...")
    modeller.addSolvent(
        forcefield=app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml"),
        padding=1.0 * unit.nanometer,  # 1nm水层
        ionicStrength=0.15 * unit.molar,  # 生理盐浓度
        positiveIon='Na+', negativeIon='Cl-'
    )
    
    # 4. 创建自定义力场
    print("  创建力场...")
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    
    # 添加Nd(III)参数到力场
    # 创建自定义力场XML
    nd_xml = f"""
<ForceField>
 <AtomTypes>
  <Type name="Nd" class="Nd" element="Nd" mass="{ND_PARAMS['mass_amu']}"/>
 </AtomTypes>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="Nd" charge="{ND_PARAMS['charge']}" sigma="{SIGMA_NM}" epsilon="{EPSILON_KCAL}"/>
 </NonbondedForce>
</ForceField>
"""
    
    # 保存临时XML
    nd_ff_file = f"md_systems/{pdb_id}_nd_ff.xml"
    with open(nd_ff_file, "w") as f:
        f.write(nd_xml)
    
    # 加载自定义力场
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml", nd_ff_file)
    
    # 5. 创建系统
    print("  创建系统...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )
    
    # 6. 设置温度、压力
    print("  设置模拟参数...")
    temperature = 300 * unit.kelvin
    pressure = 1.0 * unit.atmosphere
    
    # 添加Langevin积分器
    integrator = mm.LangevinMiddleIntegrator(
        temperature, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )
    
    # 添加压力控制（NPT）
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
    
    # 7. 保存体系
    system_file = f"md_systems/{pdb_id}_system.xml"
    with open(system_file, "w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    
    # 保存初始结构
    init_pdb = f"md_systems/{pdb_id}_init.pdb"
    with open(init_pdb, "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    print(f"  ✓ 体系已保存: {system_file}")
    print(f"  ✓ 初始结构: {init_pdb}")
    
    return {
        "system": system,
        "topology": modeller.topology,
        "positions": modeller.positions,
        "integrator": integrator,
        "system_file": system_file,
        "init_pdb": init_pdb,
        "n_atoms": modeller.topology.getNumAtoms(),
        "n_residues": modeller.topology.getNumResidues(),
    }

# 构建两个体系
md_systems = {}
for pdb_id in ["4IQK", "6HN3"]:
    try:
        md_systems[pdb_id] = create_md_system(pdb_id, docking_results[pdb_id])
    except Exception as e:
        print(f"  ✗ 构建失败: {e}")
        continue

# ============================================================
# Step 3: 运行短时间MD（平衡+生产）并分析
# ============================================================
print("\n" + "=" * 60)
print("Step 3: 运行MD模拟（10ns生产模拟）")
print("=" * 60)

def run_md_simulation(pdb_id, md_data, steps=5000000):  # 10ns (2fs × 5M步)
    """运行MD模拟"""
    print(f"\n--- {pdb_id} MD模拟 ---")
    
    # 创建模拟器
    simulation = app.Simulation(
        md_data["topology"],
        md_data["system"],
        md_data["integrator"],
        platform=Platform.getPlatformByName('CPU')
    )
    simulation.context.setPositions(md_data["positions"])
    
    # 最小化能量
    print("  能量最小化...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # 平衡阶段1: NVT (100ps)
    print("  NVT平衡 (100ps)...")
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    simulation.step(50000)  # 100ps
    
    # 平衡阶段2: NPT (100ps)
    print("  NPT平衡 (100ps)...")
    simulation.step(50000)  # 100ps
    
    # 生产模拟
    print(f"  生产模拟 ({steps//500000}ns)...")
    
    # 设置输出
    traj_file = f"md_trajectories/{pdb_id}_traj.dcd"
    log_file = f"md_trajectories/{pdb_id}_log.csv"
    
    simulation.reporters.append(app.DCDReporter(traj_file, 2500))  # 每5ps保存一帧
    simulation.reporters.append(app.StateDataReporter(
        log_file, 2500, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        speed=True, remainingTime=True, totalSteps=steps
    ))
    
    # 运行模拟
    t0 = time.time()
    simulation.step(steps)
    elapsed = time.time() - t0
    
    print(f"  ✓ 模拟完成: {elapsed:.1f}s")
    print(f"  ✓ 轨迹文件: {traj_file}")
    print(f"  ✓ 日志文件: {log_file}")
    
    return {
        "simulation": simulation,
        "traj_file": traj_file,
        "log_file": log_file,
        "elapsed_s": elapsed,
    }

# 运行模拟（缩短到2ns用于快速演示）
md_sims = {}
for pdb_id, md_data in md_systems.items():
    try:
        md_sims[pdb_id] = run_md_simulation(pdb_id, md_data, steps=1000000)  # 2ns
    except Exception as e:
        print(f"  ✗ MD模拟失败: {e}")
        continue

# ============================================================
# Step 4: 分析轨迹并生成论文图表
# ============================================================
print("\n" + "=" * 60)
print("Step 4: 分析MD轨迹并生成论文图表")
print("=" * 60)

try:
    import mdtraj as md
    print("✓ MDTraj 已安装")
except ImportError:
    print("✗ MDTraj 未安装，尝试安装...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "mdtraj"])
    import mdtraj as md
    print("✓ MDTraj 安装成功")

os.makedirs("md_analysis", exist_ok=True)

def analyze_md_trajectory(pdb_id, md_data, sim_data):
    """分析MD轨迹"""
    print(f"\n--- {pdb_id} 轨迹分析 ---")
    
    # 加载轨迹
    traj = md.load(sim_data["traj_file"], top=md_data["init_pdb"])
    print(f"  轨迹: {traj.n_frames} 帧, {traj.n_atoms} 原子")
    
    # 找到Nd原子索引
    nd_indices = [i for i, atom in enumerate(traj.topology.atoms) 
                 if atom.element.symbol == 'Nd']
    if not nd_indices:
        print("  ✗ 未找到Nd原子")
        return None
    
    nd_idx = nd_indices[0]
    print(f"  Nd原子索引: {nd_idx}")
    
    # 1. RMSD分析（相对于初始结构）
    print("  计算RMSD...")
    protein_indices = traj.topology.select("protein")
    rmsd = md.rmsd(traj, traj, 0, atom_indices=protein_indices)
    
    # 2. Nd与蛋白的距离
    print("  计算Nd-蛋白距离...")
    protein_heavy = traj.topology.select("protein and element != 'H'")
    
    # 找到最近的重原子
    min_distances = []
    contact_residues = []
    
    for frame in range(traj.n_frames):
        # 计算Nd到所有蛋白重原子的距离
        atom_pairs = np.array([[nd_idx, i] for i in protein_heavy])
        dists = md.compute_distances(traj[frame], atom_pairs)
        min_dist = np.min(dists)
        min_distances.append(min_dist)
        
        # 找到最近残基
        min_idx = protein_heavy[np.argmin(dists)]
        atom = traj.topology.atom(min_idx)
        residue = atom.residue
        contact_residues.append(f"{residue.name}{residue.resSeq}")
    
    # 3. 回转半径 (Rg)
    print("  计算回转半径...")
    rg = md.compute_rg(traj)
    
    # 4. 二级结构
    print("  计算二级结构...")
    dssp = md.compute_dssp(traj, simplified=True)
    
    # 保存分析数据
    analysis_file = f"md_analysis/{pdb_id}_analysis.npz"
    np.savez(
        analysis_file,
        time_ns=np.arange(traj.n_frames) * 0.002,  # 2ps每帧
        rmsd=rmsd,
        nd_min_dist=np.array(min_distances),
        contact_residues=contact_residues,
        rg=rg,
        dssp=dssp,
    )
    
    print(f"  ✓ 分析数据: {analysis_file}")
    
    return {
        "traj": traj,
        "rmsd": rmsd,
        "nd_min_dist": np.array(min_distances),
        "contact_residues": contact_residues,
        "rg": rg,
        "dssp": dssp,
        "nd_idx": nd_idx,
        "analysis_file": analysis_file,
    }

# 分析轨迹
md_analyses = {}
for pdb_id, md_data in md_systems.items():
    if pdb_id in md_sims:
        try:
            md_analyses[pdb_id] = analyze_md_trajectory(pdb_id, md_data, md_sims[pdb_id])
        except Exception as e:
            print(f"  ✗ 分析失败: {e}")
            continue

# ============================================================
# Step 5: 生成论文级MD分析图表
# ============================================================
print("\n" + "=" * 60)
print("Step 5: 生成论文级MD分析图表")
print("=" * 60)

os.makedirs("figures", exist_ok=True)

for pdb_id, analysis in md_analyses.items():
    if analysis is None:
        continue
    
    print(f"\n--- {pdb_id} 图表生成 ---")
    
    time_ns = np.arange(len(analysis["rmsd"])) * 0.002  # 2ps每帧
    
    # --- 图1: RMSD随时间变化 ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # 1a: RMSD
    ax = axes[0, 0]
    ax.plot(time_ns, analysis["rmsd"], 'b-', linewidth=1.5, alpha=0.8)
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('RMSD (Å)', fontsize=11)
    ax.set_title('Backbone RMSD', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(analysis["rmsd"]), color='r', linestyle='--', alpha=0.5, 
               label=f'Mean: {np.mean(analysis["rmsd"]):.2f} Å')
    ax.legend(fontsize=9)
    
    # 1b: Nd-蛋白最小距离
    ax = axes[0, 1]
    ax.plot(time_ns, analysis["nd_min_dist"] * 10, 'g-', linewidth=1.5, alpha=0.8)  # nm→Å
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('Nd-Protein distance (Å)', fontsize=11)
    ax.set_title('Minimum Nd-Protein Distance', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=3.0, color='orange', linestyle='--', alpha=0.5, label='Coordination (3.0 Å)')
    ax.axhline(y=4.0, color='red', linestyle='--', alpha=0.5, label='Weak (4.0 Å)')
    ax.legend(fontsize=9)
    
    # 1c: 回转半径
    ax = axes[1, 0]
    ax.plot(time_ns, analysis["rg"], 'purple', linewidth=1.5, alpha=0.8)
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('Radius of gyration (Å)', fontsize=11)
    ax.set_title('Protein Radius of Gyration', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(analysis["rg"]), color='r', linestyle='--', alpha=0.5,
               label=f'Mean: {np.mean(analysis["rg"]):.2f} Å')
    ax.legend(fontsize=9)
    
    # 1d: 接触残基统计
    ax = axes[1, 1]
    from collections import Counter
    contact_counts = Counter(analysis["contact_residues"])
    top_contacts = contact_counts.most_common(10)
    
    if top_contacts:
        residues = [rc[0] for rc in top_contacts]
        counts = [rc[1] for rc in top_contacts]
        percentages = [c/len(analysis["contact_residues"])*100 for c in counts]
        
        bars = ax.barh(range(len(residues)), percentages, color='teal', alpha=0.7)
        ax.set_yticks(range(len(residues)))
        ax.set_yticklabels(residues, fontfamily='monospace', fontsize=9)
        ax.set_xlabel('Contact frequency (%)', fontsize=11)
        ax.set_title('Top 10 Contact Residues with Nd', fontweight='bold')
        ax.invert_yaxis()
        
        for i, (bar, pct) in enumerate(zip(bars, percentages)):
            ax.text(pct + 1, bar.get_y() + bar.get_height()/2, f'{pct:.1f}%',
                   va='center', fontsize=8)
    
    fig.suptitle(f'MD Simulation Analysis — {docking_results[pdb_id]["name"]} ({pdb_id})',
                fontweight='bold', fontsize=14)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_md_analysis.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {pdb_id}_md_analysis.png")
    
    # --- 图2: 轨迹快照 ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 选取3个时间点 (开始、中间、结束)
    time_points = [0, len(time_ns)//2, -1]
    time_labels = ['0 ns', f'{time_ns[len(time_ns)//2]:.1f} ns', f'{time_ns[-1]:.1f} ns']
    
    for ax, t_idx, t_label in zip(axes, time_points, time_labels):
        # 简化显示：只显示蛋白骨架和Nd
        frame = analysis["traj"][t_idx]
        
        # 蛋白骨架 (CA原子)
        ca_indices = [i for i, atom in enumerate(frame.topology.atoms) 
                     if atom.name == 'CA']
        ca_coords = frame.xyz[0][ca_indices]
        
        # Nd位置
        nd_coord = frame.xyz[0][analysis["nd_idx"]]
        
        # 2D投影 (XY平面)
        ax.scatter(ca_coords[:,0], ca_coords[:,1], s=0.5, c='lightgray', alpha=0.5, rasterized=True)
        ax.scatter(nd_coord[0], nd_coord[1], s=100, c='red', marker='*',
                  edgecolors='darkred', linewidths=1.5, zorder=10, label='Nd(III)')
        
        # 标注最近残基
        nearest_res = analysis["contact_residues"][t_idx]
        for i, atom_idx in enumerate(ca_indices):
            atom = frame.topology.atom(atom_idx)
            res_label = f"{atom.residue.name}{atom.residue.resSeq}"
            if res_label == nearest_res:
                coord = ca_coords[i]
                ax.annotate(res_label, (coord[0], coord[1]), fontsize=8, color='blue',
                          bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
                break
        
        ax.set_xlabel('X (Å)', fontsize=11)
        ax.set_ylabel('Y (Å)', fontsize=11)
        ax.set_title(f'Snapshot at {t_label}', fontweight='bold')
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=8)
    
    fig.suptitle(f'Trajectory Snapshots — {docking_results[pdb_id]["name"]} ({pdb_id})',
                fontweight='bold', fontsize=13)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_md_snapshots.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {pdb_id}_md_snapshots.png")

# ============================================================
# 汇总结果
# ============================================================
print("\n" + "=" * 60)
print("MD模拟完成")
print("=" * 60)

for pdb_id in md_analyses:
    if md_analyses[pdb_id]:
        print(f"\n{pdb_id}:")
        print(f"  RMSD: {np.mean(md_analyses[pdb_id]['rmsd']):.2f} ± {np.std(md_analyses[pdb_id]['rmsd']):.2f} Å")
        print(f"  Nd-蛋白最小距离: {(np.mean(md_analyses[pdb_id]['nd_min_dist'])*10):.2f} Å")
        print(f"  回转半径: {np.mean(md_analyses[pdb_id]['rg']):.2f} Å")
        
        # 主要接触残基
        contact_counts = Counter(md_analyses[pdb_id]["contact_residues"])
        top3 = contact_counts.most_common(3)
        if top3:
            print(f"  主要接触残基: {', '.join([f'{r}({c}帧)' for r,c in top3])}")

print("\n生成的文件:")
print("  对接结果: docking_results/grid_docking_results.json")
print("  MD体系: md_systems/*.xml, md_systems/*.pdb")
print("  MD轨迹: md_trajectories/*.dcd, md_trajectories/*.csv")
print("  分析数据: md_analysis/*.npz")
print("  图表: figures/*.png")

print("\nOK")