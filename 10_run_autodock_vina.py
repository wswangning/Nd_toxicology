#!/usr/bin/env python3
"""
Run AutoDock Vina for each protein-Nd³⁺ complex.
"""
import subprocess
import os
import yaml

def prepare_receptor(pdb_file):
    # Convert to pdbqt using prepare_receptor4.py from MGLTools
    pdbqt = pdb_file.replace(".pdb", ".pdbqt")
    cmd = f"prepare_receptor4.py -r {pdb_file} -o {pdbqt}"
    subprocess.run(cmd, shell=True, check=True)
    return pdbqt

def run_vina(receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z, size_x, size_y, size_z, out_prefix):
    config_txt = f"""
receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}
size_x = {size_x}
size_y = {size_y}
size_z = {size_z}
exhaustiveness = 32
num_modes = 9
energy_range = 5
out = {out_prefix}_out.pdbqt
log = {out_prefix}_log.txt
"""
    with open(f"{out_prefix}_vina.conf", "w") as f:
        f.write(config_txt)
    cmd = f"vina --config {out_prefix}_vina.conf"
    subprocess.run(cmd, shell=True, check=True)
    return f"{out_prefix}_log.txt"

def main():
    config = yaml.safe_load(open("config.yaml"))
    ligand = "data/nd.pdbqt"
    for target in config["docking_targets"]:
        receptor_pdb = f"data/pdbs/{target['pdb']}.pdb" if "pdb" in target else f"data/pdbs/AF-{target['uniprot']}.pdb"
        receptor_pdbqt = prepare_receptor(receptor_pdb)
        run_vina(receptor_pdbqt, ligand,
                 target["center_x"], target["center_y"], target["center_z"],
                 target["size_x"], target["size_y"], target["size_z"],
                 f"outputs/docking/{target['name']}")
    print("Docking finished.")

if __name__ == "__main__":
    main()