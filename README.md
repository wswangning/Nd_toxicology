# Nd(III)分子对接与MD模拟 —— 完整工作流

## 概述

本工作流针对论文《Mechanistic Insights into Neodymium Nitrate-Induced Hepatotoxicity》
中的分子对接和分子动力学模拟部分，提供了改进的、可完全复现的计算方案。

## 工作流步骤

### 1. 环境准备
```bash
# 安装Python包
pip install rdkit biopython numpy scipy matplotlib MDAnalysis openmm pdbfixer openmmforcefields

# 安装AutoDock Vina（从官网下载）
# https://github.com/ccsb-scripps/AutoDock-Vina/releases
# 将vina.exe添加到PATH
```

### 2. 运行工作流
```bash
# 进入工作目录
cd C:\Users\wangning\AppData\Roaming\Tencent\Marvis\User\oAN1i2eqkIJ6qD9R98ixvZQsjxjI\workspace\conv_19e5eb9a626_881806416319\temp

# 运行完整工作流（已部分完成）
python complete_workflow.py

# 或手动执行各步骤
# 2.1 下载蛋白结构
python download_pdb.py

# 2.2 分子对接
vina --config docking_results/4IQK_blind_dock.conf --out docking_results/4IQK_docked.pdbqt
vina --config docking_results/4IQK_cys_focused.conf --out docking_results/4IQK_cys_docked.pdbqt

# 2.3 MD模拟（需要GPU）
# 使用OpenMM脚本，完整100ns需GPU集群
```

### 3. 轨迹分析
```bash
# 在完整MD轨迹上运行
python analyze_trajectory.py --pdb md_results/4IQK_final.pdb --traj md_results/4IQK_traj.dcd --output analysis/ --protein 4IQK
```

## 关键文件说明

### 蛋白结构
- `pdb_structures/4IQK_prepared.pdb` — NFE2L2-Keap1复合物（已准备）
- `pdb_structures/6HN3_prepared.pdb` — GPX4（已准备）

### 对接配置
- `docking_results/4IQK_blind_dock.conf` — 盲对接配置
- `docking_results/4IQK_cys_focused.conf` — Cys富集区域聚焦对接
- `docking_results/4IQK_receptor.pdbqt` — 受体PDBQT文件

### 力场参数
- `md_systems/Nd_forcefield.json` — Nd(III) LJ参数（σ, ε, q）

### 分析脚本
- `analyze_trajectory.py` — 完整的MD轨迹分析框架
- `complete_workflow.py` — 完整工作流脚本（已部分运行）

## 方法学改进要点

### 1. 对接改进
- 双策略对接：盲对接 + 金属结合位点聚焦对接
- 完整记录格点盒子参数（中心坐标 + 尺寸）
- 提高exhaustiveness（64/128 vs 论文的32）
- 明确Nd(III) +3电荷处理

### 2. MD改进
- **力场参数具体化**：σ=0.263 nm, ε=0.50 kJ/mol（论文仅模糊描述）
- **平衡时间增强**：2 ns NVT + 2 ns NPT（论文仅100 ps）
- **分析完整性**：补充Rg、氢键、能量分解、MM/GBSA等论文缺失项
- **收敛性验证**：RMSD收敛曲线 + MM/PBSA收敛分析

### 3. 可复现性
- 所有参数明文化（JSON格式）
- 配置文件可直接复用
- 文献依据完整（Shannon离子半径 + Li & Merz力场参数）

## 后续步骤建议

### 短期（验证论文结果）
1. 使用提供的对接配置运行Vina，验证论文的-8.9 kcal/mol结合能
2. 在本地GPU上运行100 ns MD（使用OpenMM脚本）
3. 运行轨迹分析，获取RMSD/RMSF/Rg等定量数据

### 中期（改进研究）
1. 使用GOLD软件进行金属配位对接（替代Vina）
2. 通过MCPB.py + QM计算推导更精确的Nd(III)参数
3. 运行QM/MM模拟，考虑极化效应

### 长期（扩展研究）
1. 扩展到全部6个靶蛋白（MT1, MT2, ACSL4, SLC7A11）
2. 进行自由能微扰（FEP）计算，获得更精确的ΔΔG
3. 结合实验验证（ITC, SPR, 晶体结构）

## 注意事项

1. **计算资源**：完整100 ns MD对每个复合物约需24-48 GPU小时
2. **软件局限**：Vina对金属离子的固有局限，建议使用GOLD
3. **力场局限**：经典LJ势无法描述极化效应，建议AMOEBA/Drude力场
4. **收敛性**：MM/PBSA结果需通过收敛性分析验证

## 联系信息

如需进一步的技术支持或合作，请提供：
- 完整的MD轨迹文件
- 对接结果
- 具体的分析需求

本工作流为论文《Mechanistic Insights into Neodymium Nitrate-Induced Hepatotoxicity》
提供了完全可复现的计算化学改进方案，显著提升了方法学的严谨性和可验证性。