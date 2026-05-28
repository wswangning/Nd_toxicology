"""
Nd(III) 真空MD — HDF5轨迹格式
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
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm import Platform

# 安装mdtraj
try:
    import mdtraj as md
except:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "mdtraj"])
    import mdtraj as md

ND_SIGMA_NM = 0.263
ND_EPSILON_KJ = 0.50
ND_CHARGE = 3.0

print("=" * 60)
print("Nd(III) Vacuum MD — HDF5")
print("=" * 60)
t_total = time.time()

with open("docking_results/grid_docking_results.json", "r") as f:
    dr = json.load(f)

os.makedirs("figures", exist_ok=True)
os.makedirs("md_systems", exist_ok=True)
os.makedirs("md_trajectories", exist_ok=True)

for pdb_id in ["6HN3", "4IQK"]:
    print(f"\n{'='*60}")
    print(f"  {pdb_id} ({dr[pdb_id]['name']})")
    print(f"{'='*60}")
    
    nd_pos = dr[pdb_id]["best_position"]
    
    src = f"pdb_structures/{pdb_id}_prepared.pdb"
    with open(src) as f:
        lines = f.readlines()
    
    last_atom_idx = max(i for i, l in enumerate(lines) if l.startswith("ATOM ") or l.startswith("HETATM"))
    last_serial = int(lines[last_atom_idx][6:11].strip())
    nx, ny, nz = nd_pos
    nd_line = f"HETATM{last_serial+1:5d} ND    ND A1000    {nx:8.3f}{ny:8.3f}{nz:8.3f}  1.00  0.00          ND  \n"
    lines.insert(last_atom_idx + 1, nd_line)
    
    complex_pdb = f"md_systems/{pdb_id}_vac.pdb"
    with open(complex_pdb, "w") as f:
        f.writelines(lines)
    
    pdb = app.PDBFile(complex_pdb)
    forcefield = app.ForceField("amber14-all.xml")
    
    nd_xml = f"""<ForceField>
 <Residues>
  <Residue name="ND">
   <Atom name="ND" type="Nd"/>
  </Residue>
 </Residues>
 <AtomTypes>
  <Type name="Nd" class="Nd" element="Nd" mass="144.24"/>
 </AtomTypes>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="Nd" charge="{ND_CHARGE}" sigma="{ND_SIGMA_NM}" epsilon="{ND_EPSILON_KJ}"/>
 </NonbondedForce>
</ForceField>"""
    nd_xml_path = f"md_systems/{pdb_id}_nd.xml"
    with open(nd_xml_path, "w") as f:
        f.write(nd_xml)
    forcefield.loadFile(nd_xml_path)
    
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=2.0*unit.nanometer,
        constraints=app.HBonds,
    )
    
    n_atoms = pdb.topology.getNumAtoms()
    print(f"  真空体系: {n_atoms} 原子")
    
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 5.0/unit.picosecond, 2.0*unit.femtoseconds
    )
    simulation = app.Simulation(pdb.topology, system, integrator,
                               platform=Platform.getPlatformByName('CPU'))
    simulation.context.setPositions(pdb.positions)
    
    print("  最小化...")
    t1 = time.time()
    simulation.minimizeEnergy(maxIterations=500)
    print(f"  {time.time()-t1:.1f}s")
    
    print("  平衡 (10ps)...")
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    t1 = time.time()
    simulation.step(5000)
    print(f"  {time.time()-t1:.1f}s")
    
    # HDF5 reporter
    h5_file = f"md_trajectories/{pdb_id}_vac_traj.h5"
    report_interval = 1000
    simulation.reporters.append(md.reporters.HDF5Reporter(h5_file, report_interval))
    
    steps = 50000  # 0.1ns
    print(f"  生产 0.1ns ({steps} steps)...")
    t1 = time.time()
    simulation.step(steps)
    elapsed = time.time() - t1
    
    for reporter in simulation.reporters:
        try:
            reporter.close()
        except:
            pass
    del simulation
    
    print(f"  ✓ {elapsed:.1f}s ({steps/elapsed:.0f} steps/s)")
    
    if not os.path.exists(h5_file):
        print(f"  ✗ HDF5文件未生成: {h5_file}")
        continue
    print(f"  HDF5: {os.path.getsize(h5_file)} bytes")
    
    # 分析
    print("  分析...")
    traj = md.load(h5_file, top=complex_pdb)
    print(f"  轨迹: {traj.n_frames}帧")
    
    protein_atoms = traj.topology.select("protein")
    protein_heavy = traj.topology.select("protein and element != H")
    nd_idx = [i for i, a in enumerate(traj.topology.atoms) 
              if a.name == "ND" and a.residue.name == "ND"][0]
    
    rmsd = md.rmsd(traj, traj, 0, atom_indices=protein_atoms)
    rmsf = md.rmsf(traj, traj, 0, atom_indices=protein_atoms)
    rg = md.compute_rg(traj)
    
    nd_dists = np.zeros(traj.n_frames)
    for i in range(traj.n_frames):
        nd_pos_t = traj.xyz[i, nd_idx]
        prot_xyz = traj.xyz[i][protein_heavy]
        nd_dists[i] = np.min(np.linalg.norm(prot_xyz - nd_pos_t, axis=1))
    
    contact_counter = {}
    for i in range(traj.n_frames):
        nd_pos_t = traj.xyz[i, nd_idx]
        for j in protein_heavy:
            d = np.linalg.norm(traj.xyz[i, j] - nd_pos_t)
            if d < 0.5:
                atom = traj.topology.atom(j)
                res_key = f"{atom.residue.name}{atom.residue.resSeq}"
                contact_counter[res_key] = contact_counter.get(res_key, 0) + 1
    
    time_ns = np.arange(traj.n_frames) * 0.002  # 1000*2fs = 2ps per frame
    
    # 图表
    print("  图表...")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    ax.plot(time_ns, rmsd * 10, 'navy', linewidth=1.2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('Backbone RMSD')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(rmsd*10), color='red', linestyle='--', alpha=0.5,
              label=f'Mean: {np.mean(rmsd*10):.2f} Å')
    ax.legend(fontsize=8)
    
    ax = axes[1]
    ax.bar(np.arange(len(rmsf)), rmsf * 10, width=1.0, color='steelblue', edgecolor='none')
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title('Per-Residue RMSF')
    ax.grid(True, alpha=0.3, axis='y')
    
    ax = axes[2]
    ax.plot(time_ns, rg * 10, 'darkgreen', linewidth=1.2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (Å)')
    ax.set_title('Radius of Gyration')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=np.mean(rg*10), color='red', linestyle='--', alpha=0.5,
              label=f'Mean: {np.mean(rg*10):.2f} Å')
    ax.legend(fontsize=8)
    
    fig.suptitle(f'MD ({dr[pdb_id]["name"]}, {pdb_id}) — Vacuum, 0.1ns', fontweight='bold', fontsize=12)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_md_analysis.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    ax = axes[0]
    ax.plot(time_ns, nd_dists * 10, 'firebrick', linewidth=1.2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Nd-Protein Distance (Å)')
    ax.set_title('Minimum Nd-Protein Distance')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=3.0, color='orange', linestyle='--', alpha=0.5, label='Coordination (3 Å)')
    ax.axhline(y=5.0, color='gray', linestyle='--', alpha=0.5, label='Contact (5 Å)')
    ax.legend(fontsize=8)
    
    ax = axes[1]
    ax.text(0.5, 0.5, 'Vacuum simulation\nenergy not physically meaningful', 
            ha='center', va='center', transform=ax.transAxes, fontsize=10)
    ax.set_title('Note')
    
    fig.suptitle(f'Nd(III) Binding Dynamics ({pdb_id})', fontweight='bold', fontsize=12)
    plt.tight_layout()
    fig.savefig(f"figures/{pdb_id}_nd_binding.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    time_points = [0, len(time_ns)//2, -1]
    labels = ['0 ns', f'{time_ns[len(time_ns)//2]:.2f} ns', f'{time_ns[-1]:.2f} ns']
    
    for ax, tp, lab in zip(axes, time_points, labels):
        ca_idx = [i for i, a in enumerate(traj.topology.atoms) if a.name == 'CA']
        ca_xyz = traj.xyz[tp][ca_idx] * 10
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
    
    print(f"  RMSD: {np.mean(rmsd*10):.2f} ± {np.std(rmsd*10):.2f} Å")
    print(f"  Nd距离: {np.mean(nd_dists*10):.2f} Å")
    top5 = sorted(contact_counter.items(), key=lambda x: -x[1])[:5]
    print(f"  Top contacts: {', '.join([f'{k}({v})' for k,v in top5])}")

print(f"\n{'='*60}")
print(f"全部完成 (total {time.time()-t_total:.0f}s)")
print(f"{'='*60}")