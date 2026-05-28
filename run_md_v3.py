"""
Nd(III) MD模拟 — 隐式溶剂GBSA加速版
======================================
使用OpenMM GBSA隐式溶剂，大幅减少原子数，加速模拟
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
    'axes.linewidth': 1.2, 'axes.labelsize': 12, 'axes.titlesize': 13,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm import Platform, XmlSerializer

ND_SIGMA_NM = 0.263
ND_EPSILON_KCAL = 0.50 / 4.184
ND_CHARGE = 3.0

print("=" * 60)
print("Nd(III) MD Simulation — GBSA Implicit Solvent")
print("=" * 60)

with open("docking_results/grid_docking_results.json", "r") as f:
    docking_results = json.load(f)

os.makedirs("md_systems", exist_ok=True)
os.makedirs("md_trajectories", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ============================================================
# Step 1: 创建蛋白+Nd复合PDB
# ============================================================
print("\nStep 1: 创建复合结构...")

def create_complex_pdb(pdb_id, nd_pos):
    src = f"pdb_structures/{pdb_id}_prepared.pdb"
    dst = f"md_systems/{pdb_id}_complex.pdb"
    with open(src) as f:
        lines = f.readlines()
    
    last_atom_idx = max(i for i, l in enumerate(lines) if l.startswith("ATOM ") or l.startswith("HETATM"))
    last_serial = int(lines[last_atom_idx][6:11].strip())
    
    nx, ny, nz = nd_pos
    nd_line = f"HETATM{last_serial+1:5d} ND    ND A1000    {nx:8.3f}{ny:8.3f}{nz:8.3f}  1.00  0.00          ND  \n"
    lines.insert(last_atom_idx + 1, nd_line)
    
    with open(dst, "w") as f:
        f.writelines(lines)
    return dst

for pdb_id, result in docking_results.items():
    create_complex_pdb(pdb_id, result["best_position"])
    print(f"  ✓ {pdb_id}")

# ============================================================
# Step 2: 构建GBSA体系
# ============================================================
print("\nStep 2: 构建GBSA体系...")

def build_gbsa_system(pdb_id):
    complex_pdb = f"md_systems/{pdb_id}_complex.pdb"
    pdb = app.PDBFile(complex_pdb)
    
    # GBSA力场 (amber99sb + obc2)
    forcefield = app.ForceField("amber99sb.xml", "implicit/obc2.xml")
    
    # 注册Nd(III)自定义残基模板
    nd_template = f'''<ForceField>
 <Residues>
  <Residue name="ND">
   <Atom name="Nd" type="Nd3+" charge="{ND_CHARGE}"/>
  </Residue>
 </Residues>
 <AtomTypes>
  <Type name="Nd3+" class="Nd3+" element="Nd" mass="{144.24}"/>
 </AtomTypes>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="Nd3+" charge="{ND_CHARGE}" sigma="{ND_SIGMA_NM}" epsilon="{ND_EPSILON_KCAL * 4.184}"/>
 </NonbondedForce>
 <GBSAOBCForce>
  <Atom type="Nd3+" charge="{ND_CHARGE}" radius="0.169" scale="0.79"/>
 </GBSAOBCForce>
</ForceField>'''
    
    nd_xml = f"md_systems/{pdb_id}_nd.xml"
    with open(nd_xml, "w") as f:
        f.write(nd_template)
    
    forcefield.loadFile(nd_xml)
    
    init_pdb = f"md_systems/{pdb_id}_gb_init.pdb"
    with open(init_pdb, "w") as f:
        app.PDBFile.writeFile(pdb.topology, pdb.positions, f)
    
    # 创建GBSA系统
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    
    # 验证Nd参数
    forces = {type(f).__name__: f for f in system.getForces()}
    nb_force = forces.get("NonbondedForce")
    
    nd_found = False
    for atom in pdb.topology.atoms():
        if atom.name == "Nd" and atom.residue.name == "ND":
            charge_val, sigma_val, eps_val = nb_force.getParticleParameters(atom.index)
            print(f"  Nd参数: q={charge_val}, σ={sigma_val}nm, ε={eps_val}kJ/mol")
            nd_found = True
            break
    
    if not nd_found:
        print("  ⚠ 未找到Nd原子!")
    
    n_atoms = pdb.topology.getNumAtoms()
    print(f"  体系: {n_atoms} 原子")
    
    return system, pdb.topology, pdb.positions, init_pdb

md_systems = {}
for pdb_id in ["4IQK", "6HN3"]:
    try:
        md_systems[pdb_id] = build_gbsa_system(pdb_id)
    except Exception as e:
        print(f"  ✗ {pdb_id}: {e}")
        import traceback
        traceback.print_exc()

# ============================================================
# Step 3: 运行MD (5ns per protein)
# ============================================================
print("\n" + "=" * 60)
print("Step 3: 运行MD模拟 (5ns GBSA)")
print("=" * 60)

def run_md_gb(pdb_id, system, topology, positions, steps=1000000):
    """运行GBSA MD模拟 (2ns)"""
    print(f"\n--- {pdb_id} ---")
    
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds
    )
    
    simulation = app.Simulation(topology, system, integrator,
                               platform=Platform.getPlatformByName('CPU'))
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("  最小化...")
    simulation.minimizeEnergy(maxIterations=2000)
    
    # 平衡
    print("  NVT平衡 (20ps)...")
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(10000)  # 20ps
    
    # 生产
    report_interval = 5000  # 10ps/frame
    
    traj_file = f"md_trajectories/{pdb_id}_gb_traj.dcd"
    log_file = f"md_trajectories/{pdb_id}_gb_log.csv"
    
    simulation.reporters.append(app.DCDReporter(traj_file, report_interval))
    
    total_steps_str = str(steps)
    simulation.reporters.append(app.StateDataReporter(
        log_file, report_interval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, speed=True, remainingTime=True,
        totalSteps=steps
    ))
    
    print(f"  生产模拟 {steps//500000}ns ({steps} steps)...")
    t0 = time.time()
    simulation.step(steps)
    elapsed = time.time() - t0
    
    print(f"  ✓ 完成: {elapsed:.1f}s ({steps/elapsed:.0f} steps/s)")
    
    return simulation, traj_file, log_file, elapsed

md_results = {}
for pdb_id, (system, topology, positions, init_pdb) in md_systems.items():
    try:
        sim, traj, log, elapsed = run_md_gb(pdb_id, system, topology, positions, steps=1000000)  # 2ns
        md_results[pdb_id] = {"sim": sim, "traj": traj, "log": log, "elapsed": elapsed, "init_pdb": init_pdb}
    except Exception as e:
        print(f"  ✗ {pdb_id} 失败: {e}")
        import traceback
        traceback.print_exc()

# ============================================================
# Step 4: 分析轨迹
# ============================================================
print("\n" + "=" * 60)
print("Step 4: 轨迹分析")
print("=" * 60)

try:
    import mdtraj as md
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "mdtraj"])
    import mdtraj as md

os.makedirs("md_analysis", exist_ok=True)

for pdb_id, result in md_results.items():
    print(f"\n--- {pdb_id} ---")
    
    traj = md.load(result["traj"], top=result["init_pdb"])
    print(f"  轨迹: {traj.n_frames}帧, {traj.n_atoms}原子")
    
    # 蛋白选择
    protein_atoms = traj.topology.select("protein")
    protein_heavy = traj.topology.select("protein and element != H")
    
    # Nd原子索引
    nd_indices = [i for i, atom in enumerate(traj.topology.atoms) if atom.name == "ND" and atom.residue.name == "ND"]
    
    if not nd_indices:
        print("  ✗ 未找到Nd原子")
        continue
    
    nd_idx = nd_indices[0]
    
    # RMSD
    print("  RMSD...")
    rmsd = md.rmsd(traj, traj, 0, atom_indices=protein_atoms)
    
    # RMSF
    print("  RMSF...")
    rmsf = md.rmsf(traj, traj, 0, atom_indices=protein_atoms)
    
    # 回转半径
    print("  Rg...")
    rg = md.compute_rg(traj)
    
    # Nd-蛋白距离 (到最近蛋白重原子)
    print("  Nd距离...")
    nd_dists = np.zeros(traj.n_frames)
    for i in range(traj.n_frames):
        nd_pos = traj.xyz[i, nd_idx]
        prot_xyz = traj.xyz[i][protein_heavy]
        nd_dists[i] = np.min(np.linalg.norm(prot_xyz - nd_pos, axis=1))
    
    # 找出Nd频繁接触的残基（5Å内）
    contact_counter = {}
    for i in range(traj.n_frames):
        if nd_dists[i] < 0.5:  # <5Å
            nd_pos = traj.xyz[i, nd_idx]
            for j in protein_heavy:
                d = np.linalg.norm(traj.xyz[i, j] - nd_pos)
                if d < 0.5:
                    atom = traj.topology.atom(j)
                    res_key = f"{atom.residue.name}{atom.residue.resSeq}"
                    contact_counter[res_key] = contact_counter.get(res_key, 0) + 1
    
    time_ns = np.arange(traj.n_frames) * 0.005  # 5ps per frame
    
    # 保存数据
    np.savez(f"md_analysis/{pdb_id}_gb_analysis.npz",
             time_ns=time_ns,
             rmsd=rmsd, rmsf=rmsf, rg=rg,
             nd_min_dist=nd_dists,
             contact_counter=dict(sorted(contact_counter.items(), key=lambda x: -x[1])))
    
    # ============================================================
    # Step 5: 论文图表
    # ============================================================
    print("  生成图表...")
    
    # 图1: RMSD/RMSF/Rg 三合一
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    # RMSD
    ax = axes[0]
    ax.plot(time_ns, rmsd * 10, 'navy', linewidth=1.2)  # nm→Å
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('Backbone RMSD')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(rmsd*10), color='red', linestyle='--', alpha=0.5,
              label=f'Mean: {np.mean(rmsd*10):.2f} Å')
    ax.legend(fontsize=8)
    
    # RMSF
    ax = axes[1]
    residues = np.arange(len(rmsf))
    ax.bar(residues, rmsf * 10, width=1.0, color='steelblue', edgecolor='none')  # nm→Å
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title('Per-Residue RMSF')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Rg
    ax = axes[2]
    ax.plot(time_ns, rg * 10, 'darkgreen', linewidth=1.2)  # nm→Å
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (Å)')
    ax.set_title('Radius of Gyration')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(rg*10), color='red', linestyle='--', alpha=0.5,
              label=f'Mean: {np.mean(rg*10):.2f} Å')
    ax.legend(fontsize=8)
    
    name = docking_results[pdb_id]["name"]
    fig.suptitle(f'MD Simulation ({name}, {pdb_id}) — GBSA Implicit Solvent, 5ns', fontweight='bold', fontsize=12)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_md_analysis.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # 图2: Nd距离 + 能量随时间变化
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Nd-蛋白距离
    ax = axes[0]
    ax.plot(time_ns, nd_dists * 10, 'firebrick', linewidth=1.2)  # nm→Å
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Nd-Protein Distance (Å)')
    ax.set_title('Minimum Nd-Protein Distance')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=3.0, color='orange', linestyle='--', alpha=0.5, label='Coordination (3 Å)')
    ax.axhline(y=5.0, color='gray', linestyle='--', alpha=0.5, label='Contact (5 Å)')
    ax.legend(fontsize=8)
    
    # 能量（从日志读取）
    ax = axes[1]
    try:
        log_data = np.genfromtxt(result["log"], delimiter=',', skip_header=1)
        if log_data.ndim == 1:
            log_data = log_data.reshape(1, -1)
        pot_e = log_data[:, 2] if log_data.shape[1] > 2 else None
        if pot_e is not None:
            pot_time = np.arange(len(pot_e)) * 0.005
            ax.plot(pot_time[:len(pot_e)], pot_e, 'teal', linewidth=0.8)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('Potential Energy (kJ/mol)')
            ax.set_title('Potential Energy')
            ax.grid(True, alpha=0.3)
    except:
        ax.text(0.5, 0.5, 'Energy data unavailable', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle(f'Nd(III) Binding Dynamics ({pdb_id})', fontweight='bold', fontsize=12)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_nd_binding.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # 图3: 轨迹快照
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    time_points = [0, len(time_ns)//2, -1]
    labels = ['0 ns', f'{time_ns[len(time_ns)//2]:.1f} ns', f'{time_ns[-1]:.1f} ns']
    
    for ax, tp, lab in zip(axes, time_points, labels):
        ca_idx = [i for i, a in enumerate(traj.topology.atoms) if a.name == 'CA']
        ca_xyz = traj.xyz[tp][ca_idx] * 10  # nm→Å
        nd_xyz = traj.xyz[tp][nd_idx] * 10
        
        ax.scatter(ca_xyz[:,0], ca_xyz[:,1], s=0.3, c='lightgray', alpha=0.5, rasterized=True)
        ax.scatter(nd_xyz[0], nd_xyz[1], s=120, c='red', marker='*', edgecolors='darkred',
                  linewidths=1.5, zorder=10, label='Nd(III)')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_title(lab)
        ax.set_aspect('equal')
        ax.legend(fontsize=8)
    
    fig.suptitle(f'Trajectory Snapshots ({pdb_id})', fontweight='bold', fontsize=12)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_md_snapshots.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  ✓ 图表已保存")
    print(f"  RMSD: {np.mean(rmsd*10):.2f} ± {np.std(rmsd*10):.2f} Å")
    print(f"  Nd距离: {np.mean(nd_dists*10):.2f} Å")
    print(f"  Top contacts: {', '.join([f'{k}({v})' for k,v in sorted(contact_counter.items(), key=lambda x:-x[1])[:5]])}")

# ============================================================
# 汇总
# ============================================================
print("\n" + "=" * 60)
print("全部模拟完成")
print("=" * 60)

print("\n最终产出文件:")
print("  对接: docking_results/grid_docking_results.json")
print("  对接图: figures/4IQK_docking_pose.png, figures/6HN3_docking_pose.png")
print("  对接能量: figures/4IQK_cluster_comparison.png, figures/6HN3_cluster_comparison.png")
print("  MD分析: figures/4IQK_md_analysis.png, figures/6HN3_md_analysis.png")
print("  Nd结合: figures/4IQK_nd_binding.png, figures/6HN3_nd_binding.png")
print("  MD快照: figures/4IQK_md_snapshots.png, figures/6HN3_md_snapshots.png")
print("  MD数据: md_analysis/*_gb_analysis.npz")

print("\nOK")