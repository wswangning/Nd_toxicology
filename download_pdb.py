import requests
import os

# 创建目录
os.makedirs("pdb_structures", exist_ok=True)
os.makedirs("ligands", exist_ok=True)

# 论文中提到的PDB ID
pdb_ids = ["4IQK", "6HN3"]  # NFE2L2-Keap1 和 GPX4
# MT1, MT2, ACSL4, SLC7A11 (AlphaFold2) 稍后处理

for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(f"pdb_structures/{pdb_id}.pdb", "w") as f:
            f.write(response.text)
        print(f"下载成功: {pdb_id}")
    except Exception as e:
        print(f"下载失败 {pdb_id}: {e}")

# 创建Nd(III)离子结构文件 (PDB格式)
nd_ion_pdb = """HETATM    1  ND  ND     1       0.000   0.000   0.000  1.00  0.00          ND
END"""

with open("ligands/Nd_ion.pdb", "w") as f:
    f.write(nd_ion_pdb)

print("Nd(III)离子结构文件已创建")
print("\n文件结构:")
for root, dirs, files in os.walk("."):
    level = root.replace(".", "").count(os.sep)
    indent = " " * 2 * level
    print(f"{indent}{os.path.basename(root)}/")
    subindent = " " * 2 * (level + 1)
    for file in files:
        if file.endswith((".py", ".pdb", ".txt")):
            print(f"{subindent}{file}")