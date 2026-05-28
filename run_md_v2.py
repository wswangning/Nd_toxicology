"""
Nd(III) MD模拟 — 修正版
策略：将Nd离子作为HETATM写入PDB文件，统一加载
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
    'axes.linewidth': 1.2, 'axes.labelsize': 12, 'axes.titlesize': 13,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm import Platform, XmlSerializer

# ============================================================
# Nd(III) 力场参数
# ============================================================
ND_SIGMA_NM = 0.263
ND_EPSILON_KCAL = 0.50 / 4.184
ND_CHARGE = 3.0
ND_MASS = 144.24

print("=" * 60)
print("Step 1: 准备含Nd的复合结构")
print("=" * 60)

# 读取对接结果
with open("docking_results/grid_docking_results.json", "r") as f:
    docking_results = json.load(f)

os.makedirs("md_systems", exist_ok=True)
os.makedirs("md_trajectories", exist_ok=True)

def create_complex_pdb(pdb_id, nd_pos):
    """创建蛋白+Nd的复合PDB"""
    src = f"pdb_structures/{pdb_id}_prepared.pdb"
    dst = f"md_systems/{pdb_id}_complex.pdb"
    
    with open(src, "r") as fin:
        lines = fin.readlines()
    
    # 找到最后一个ATOM行
    last_atom_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("ATOM ") or line.startswith("HETATM"):
            last_atom_idx = i
    
    if last_atom_idx == -1:
        raise ValueError("No ATOM record found")
    
    # 解析最后一个ATOM的序列号
    last_line = lines[last_atom_idx]
    last_serial = int(last_line[6:11].strip())
    
    # 添加Nd HETATM行
    nd_x, nd_y, nd_z = nd_pos
    nd_line = (
        f"HETATM{last_serial+1:5d} ND    ND A1000    "
        f"{nd_x:8.3f}{nd_y:8.3f}{nd_z:8.3f}"
        f"  1.00  0.00          ND  \n"
    )
    
    lines.insert(last_atom_idx + 1, nd_line)
    
    # 更新TER和END
    with open(dst, "w") as fout:
        fout.writelines(lines)
    
    return dst

# 创建复合结构
for pdb_id, result in docking_results.items():
    nd_pos = result["best_position"]
    complex_pdb = create_complex_pdb(pdb_id, nd_pos)
    print(f"  ✓ {complex_pdb}")

# ============================================================
# Step 2: 构建MD体系
# ============================================================
print("\n" + "=" * 60)
print("Step 2: 构建MD体系")
print("=" * 60)

def build_system(pdb_id):
    """构建MD体系"""
    complex_pdb = f"md_systems/{pdb_id}_complex.pdb"
    
    print(f"\n{pdb_id}: 加载复合结构...")
    pdb = app.PDBFile(complex_pdb)
    
    print(f"  拓扑: {pdb.topology.getNumAtoms()} 原子, {pdb.topology.getNumResidues()} 残基")
    
    # 创建力场
    # 先做一次标准力场尝试，如果要添加Nd自定义参数需要额外处理
    print("  创建AMBER力场...")
    
    try:
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    except:
        forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3pfb.xml")
    
    # 创建Modeller
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # 添加水盒子
    print("  添加水盒子...")
    modeller.addSolvent(
        forcefield,
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        positiveIon='Na+', negativeIon='Cl-'
    )
    
    print(f"  溶剂化后: {modeller.topology.getNumAtoms()} 原子")
    
    # 保存初始结构
    init_pdb = f"md_systems/{pdb_id}_solvated.pdb"
    with open(init_pdb, "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    # 创建系统
    print("  创建系统...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    
    # 添加自定义Nd参数
    # 获取NonbondedForce
    forces = {type(f).__name__: f for f in system.getForces()}
    nb_force = forces.get("NonbondedForce")
    
    if nb_force is not None:
        # 找到Nd原子对应的粒子索引
        for atom in modeller.topology.atoms():
            if atom.name == "ND" and atom.residue.name == "ND":
                nd_index = atom.index
                # 打印当前参数
                charge_val, sigma_val, epsilon_val = nb_force.getParticleParameters(nd_index)
                print(f"  Nd当前参数: charge={charge_val}, sigma={sigma_val}nm, epsilon={epsilon_val}kJ/mol")
                # 修改为正确的Nd参数
                nb_force.setParticleParameters(nd_index, 
                    ND_CHARGE * unit.elementary_charge,
                    ND_SIGMA_NM * unit.nanometer,
                    ND_EPSILON_KCAL * unit.kilocalorie_per_mole / unit.kilojoule_per_mole * unit.kilojoule_per_mole
                )
                print(f"  Nd修改后: charge={ND_CHARGE}e, sigma={ND_SIGMA_NM}nm, epsilon={ND_EPSILON_KCAL:.4f}kcal/mol")
                break
    
    # 保存系统
    system_file = f"md_systems/{pdb_id}_system.xml"
    with open(system_file, "w") as f:
        f.write(XmlSerializer.serialize(system))
    
    print(f"  ✓ 体系已保存")
    
    return {
        "system": system,
        "topology": modeller.topology,
        "positions": modeller.positions,
        "system_file": system_file,
        "init_pdb": init_pdb,
        "n_atoms": modeller.topology.getNumAtoms(),
    }

md_systems = {}
for pdb_id in ["4IQK", "6HN3"]:
    try:
        md_systems[pdb_id] = build_system(pdb_id)
    except Exception as e:
        print(f"  ✗ 构建失败: {e}")
        import traceback
        traceback.print_exc()

# ============================================================
# Step 3: 运行MD模拟
# ============================================================
print("\n" + "=" * 60)
print("Step 3: 运行MD模拟 (2ns)")
print("=" * 60)

def run_md(pdb_id, md_data, steps=1000000):
    """运行MD模拟"""
    print(f"\n--- {pdb_id} ---")
    
    simulation = app.Simulation(
        md_data["topology"], md_data["system"],
        mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds),
        platform=Platform.getPlatformByName('CPU')
    )
    simulation.context.setPositions(md_data["positions"])
    
    # 能量最小化
    print("  能量最小化...")
    simulation.minimizeEnergy(maxIterations=2000)
    
    # 平衡
    print("  NVT平衡 (50ps)...")
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    simulation.step(25000)  # 50ps
    
    print("  NPT平衡 (50ps)...")
    simulation.step(25000)
    
    # 生产
    print(f"  生产模拟 ({steps//500000}ns)...")
    traj_file = f"md_trajectories/{pdb_id}_traj.dcd"
    log_file = f"md_trajectories/{pdb_id}_log.csv"
    
    simulation.reporters.append(app.DCDReporter(traj_file, 5000))  # 每10ps
    simulation.reporters.append(app.StateDataReporter(
        log_file, 5000, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        speed=True, remainingTime=True, totalSteps=steps
    ))
    
    t0 = time.time()
    simulation.step(steps)
    elapsed = time.time() - t0
    
    print(f"  ✓ 完成: {elapsed:.1f}s")
    return {
        "simulation": simulation, "traj_file": traj_file,
        "log_file": log_file, "elapsed": elapsed
    }

md_results = {}
for pdb_id, md_data in md_systems.items():
    try:
        md_results[pdb_id] = run_md(pdb_id, md_data, steps=1000000)
    except Exception as e:
        print(f"  ✗ 失败: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "=" * 60)
print("MD模拟完成")
print("=" * 60)

for pdb_id, r in md_results.items():
    print(f"\n{pdb_id}: {r['elapsed']:.1f}s")
    print(f"  轨迹: {r['traj_file']}")
    print(f"  日志: {r['log_file']}")

print("\nOK")