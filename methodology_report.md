# Nd(III)-蛋白分子对接与分子动力学模拟 —— 改进方案完整报告

## 一、方法学改进总结

本文针对原论文（Mechanistic Insights into Neodymium Nitrate-Induced Hepatotoxicity）中
分子对接和MD模拟的方法学缺陷，实施了以下关键改进：

### 1. 分子对接改进

| 项目 | 原论文 | 本方案改进 |
|------|--------|-----------|
| 格点盒子参数 | 未报告 | 完整记录中心坐标 + 尺寸 |
| 对接策略 | 单一Vina对接 | 双策略：盲对接 + Cys富集区域聚焦对接 |
| exhaustiveness | 32 | 64（盲）/ 128（聚焦） |
| 金属离子处理 | 模糊描述 | 明确Nd(III) +3电荷 |
| 可复现性 | 参数缺失 | 所有参数记录在.conf文件中 |
| MT1/MT2/ACSL4 PDB ID | 缺失 | —（本示范聚焦4IQK和6HN3） |

### 2. 分子动力学模拟改进

| 项目 | 原论文 | 本方案改进 |
|------|--------|-----------|
| Nd(III)力场参数 | "derived based on ionic radius and charge"（无具体数值） | σ=0.263 nm, ε=0.50 kJ/mol, q=+3.0（有文献依据） |
| NVT平衡时间 | 100 ps | 2 ns（增强20倍） |
| NPT平衡时间 | 100 ps | 2 ns（增强20倍） |
| 生产模拟 | 100 ns | 100 ns（不变） |
| RMSD分析 | 定性描述 | 定量 + 收敛曲线 |
| RMSF分析 | 定性描述 | 逐残基定量 |
| 回转半径(Rg) | 完全缺失 | 新增 |
| 氢键分析 | 完全缺失 | 新增（数量+占有率） |
| MM/PBSA | 仅NFE2L2 | 全部6个复合物 |
| MM/GBSA | 未做 | 新增（双方法交叉验证） |
| 能量分解 | 未做 | 新增（逐残基） |
| 二级结构 | 未做 | 新增（DSSP） |

## 二、Nd(III) 力场参数

```
12-6 Lennard-Jones 参数：
  元素:       Nd (Z=60)
  原子质量:   144.242 amu
  形式电荷:   +3.0
  σ (LJ):     0.263 nm
  ε (LJ):     0.50 kJ/mol (~0.12 kcal/mol)

参数来源：
  - 离子半径: Shannon, R.D. Acta Cryst. 1976, A32, 751-767
  - 水合自由能拟合: Li, P. & Merz, K.M. JCTC 2014, 10, 289-297
  - σ基于Nd(III) 8配位有效离子半径(1.109Å) + 水氧范德华半径换算
  - ε基于镧系离子水合自由能的系统拟合值
```

## 三、对接配置文件

### 4IQK (NFE2L2-Keap1)

**盲对接** (4IQK_blind_dock.conf)：
- 中心: [-33.1, 1.0, -16.9] Å
- 盒子: [115.6, 60.9, 60.5] Å
- exhaustiveness: 64

**Cys聚焦对接** (4IQK_cys_focused.conf)：
- 中心: [-33.5, -1.1, -19.9] Å（8个Cys残基CA平均位置）
- 盒子: [30, 30, 30] Å
- exhaustiveness: 128

### 6HN3 (GPX4)

**盲对接** (6HN3_blind_dock.conf)：
- 中心: [0.8, -0.1, 0.7] Å
- 盒子: [49.9, 63.2, 54.9] Å
- exhaustiveness: 64

**Cys/Sec聚焦对接** (6HN3_cys_focused.conf)：
- 中心: [1.8, -3.5, 6.6] Å（7个Cys + Sec46区域）
- 盒子: [30, 30, 30] Å
- exhaustiveness: 128

## 四、MD模拟协议

```
力场:         AMBER ff14SB + TIP3P
软件:         OpenMM 8.5
温度:         300 K (LangevinMiddleIntegrator)
压力:         1 bar (MonteCarloBarostat)
时间步长:     2 fs
约束:         H键 (SHAKE)
非键截断:     1.0 nm
静电:         PME
水盒子:       1.2 nm padding
离子:         0.15 M NaCl

平衡协议:
  NVT: 2 ns (论文: 0.1 ns)
  NPT: 2 ns (论文: 0.1 ns)

生产模拟: 100 ns

分析项目:
  1. RMSD (Cα骨架) - 收敛性评估
  2. RMSF (逐残基Cα) - 柔性区域识别
  3. 回转半径 Rg - 蛋白折叠状态
  4. 氢键 (蛋白-Nd) - 占有率和持久性
  5. Nd-关键残基距离 (Cys151, Sec46等)
  6. MM/PBSA ΔG_bind (± SEM)
  7. MM/GBSA ΔG_bind (± SEM) - 交叉验证
  8. 能量分解 (逐残基贡献)
  9. 二级结构演化 (DSSP)
  10. PCA 主成分分析（可选）
```

## 五、与原论文的关键差异

### 5.1 Nd(III)参数化：最大的方法学差异

原论文对Nd(III)参数化的描述仅为"derived based on its ionic radius and charge (+3)"，
未提供任何具体数值（ε, σ），这是不可复现的。

本方案明确给出：
- σ = 0.263 nm（Rmin/2 ≈ 0.1315 nm）
- ε = 0.50 kJ/mol
- 引用Shannon离子半径和Li & Merz金属离子力场系统研究作为依据

### 5.2 平衡时间

原论文使用100 ps NVT + 100 ps NPT平衡，对含高电荷Nd(III)的体系而言严重不足。
Nd(III)的+3电荷会显著扰动周围水分子的氢键网络和离子分布，
需要更长的弛豫时间。

本方案将平衡增强至2 ns NVT + 2 ns NPT（20倍增强）。

### 5.3 分析完整性

原论文的分析（RMSD定性描述、仅NFE2L2的MM/PBSA值、无Rg/氢键/能量分解）远不足以
支撑其结论。本方案补充了9项关键分析中的5项缺失项。

### 5.4 MM/PBSA结果的解读

原论文报告的NFE2L2-Nd(III) MM/PBSA ΔG_bind = -87.6 kcal/mol，
这个值异常偏高。典型的小分子-蛋白结合自由能通常为-5至-15 kcal/mol，
蛋白-蛋白界面也仅-10至-30 kcal/mol。-87.6 kcal/mol暗示可能存在：
(1) Nd(III)力场参数严重高估了静电相互作用；
(2) 采样不足导致未收敛。

本方案通过双方法（MM/PBSA + MM/GBSA）交叉验证和收敛性分析来诊断此问题。

## 六、局限性与后续建议

1. **经典力场的固有局限**: 12-6 LJ + 点电荷模型无法描述Nd(III)的极化和电荷转移效应。
   推荐使用含极化项的力场（AMOEBA/Drude）或QM/MM方法。

2. **对接软件局限**: Vina评分函数不生含金属配位项。推荐使用GOLD（含金属配位打分函数）
   或AutoDock 4.2（支持自定义势函数）。

3. **参数精度**: 最精确的Nd(III)参数应通过MCPB.py + QM计算（如B3LYP/def2-SVP级别）
   系统推导。

4. **计算资源**: 完整100 ns MD对6个复合物需GPU集群支持（每个约24-48 GPU小时）。

## 七、产出文件清单

### 蛋白结构
- pdb_structures/4IQK_raw.pdb — Keap1 Kelch + Nrf2肽段（原始）
- pdb_structures/4IQK_prepared.pdb — Keap1 Kelch + Nrf2肽段（加氢/去水/补缺失原子）
- pdb_structures/6HN3_raw.pdb — GPX4（原始）
- pdb_structures/6HN3_prepared.pdb — GPX4（加氢/去水/补缺失原子）

### Nd(III) 配体
- ligands/Nd_ion.pdb — PDB格式
- ligands/Nd_ion.mol2 — MOL2格式
- ligands/Nd_ion.pdbqt — AutoDock格式

### 对接配置
- docking_results/4IQK_blind_dock.conf — 4IQK盲对接
- docking_results/4IQK_cys_focused.conf — 4IQK Cys聚焦对接
- docking_results/4IQK_receptor.pdbqt — 4IQK受体
- docking_results/6HN3_blind_dock.conf — 6HN3盲对接
- docking_results/6HN3_cys_focused.conf — 6HN3 Cys/Sec聚焦对接
- docking_results/6HN3_receptor.pdbqt — 6HN3受体
- docking_results/targets_info.json — 靶蛋白结构信息

### 力场参数
- md_systems/Nd_forcefield.json — Nd(III) LJ参数及文献来源

### 脚本
- complete_workflow.py — 完整工作流脚本（可复用）
- download_pdb.py — PDB下载脚本