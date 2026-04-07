#!/usr/bin/env python3
"""
Download protein structures (PDB or AlphaFold) and prepare ligand (Nd³⁺) files.
"""
import os
import requests
import yaml
from Bio.PDB import PDBList

def download_pdb(pdb_id, outdir="data/pdbs"):
    os.makedirs(outdir, exist_ok=True)
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir=outdir, file_format="pdb")
    return os.path.join(outdir, f"{pdb_id}.pdb")

def download_alphafold(uniprot_id, outdir="data/pdbs"):
    # AlphaFold DB API
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    outpath = os.path.join(outdir, f"AF-{uniprot_id}.pdb")
    with open(outpath, "w") as f:
        f.write(response.text)
    return outpath

def prepare_nd_parameter():
    # Write a dummy PDBQT for Nd³⁺ (charge +3)
    nd_content = """REMARK  Neodymium ion
ATOM      1  ND   ND     1     0.000   0.000   0.000  1.00  0.00           ND3+
END
"""
    with open("data/nd.pdb", "w") as f:
        f.write(nd_content)
    # Use AutoDock Tools to convert to PDBQT (here we just copy)
    os.system("obabel data/nd.pdb -O data/nd.pdbqt -p 7.4")
    print("Nd³⁺ parameter prepared.")

def main():
    config = yaml.safe_load(open("config.yaml"))
    targets = config["docking_targets"]  # list of dicts: {pdb_id or uniprot, chain}
    for t in targets:
        if "pdb" in t:
            download_pdb(t["pdb"])
        elif "uniprot" in t:
            download_alphafold(t["uniprot"])
    prepare_nd_parameter()

if __name__ == "__main__":
    main()