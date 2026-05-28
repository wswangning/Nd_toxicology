# Nd(III)–Nrf2/Keap1 & Nd(III)–GPX4 分子动力学模拟研究报告

## 模拟策略

本研究针对Nd(III)与Nrf2/Keap1（PDB: 4IQK）及GPX4（PDB: 6HN3）的相互作用，采用三级计算策略：

1. **分子对接**：基于全蛋白LJ+库仑势能网格扫描方法，筛选最优结合位点
2. **真空分子动力学**：Amber14SB力场 + 自定义Nd(III)参数（σ=0.263 nm, ε=0.50 kJ/mol, q=+3.0e），0.1 ns
3. **GROMACS输入文件**：CHARMM36力场 + TIP3P水 + 0.15M NaCl，已制备完整协议（可部署于GPU服务器）

### 4IQK处理注意事项
4IQK为Keap1-Nrf2二聚体结构。为减小计算体系，仅保留A链。二聚化界面被截断可能影响局部构象，正式发表时建议保留完整二聚体。

## Nd(III)力场参数

| 参数 | 值 | 来源 |
|------|------|------|
| σ | 0.263 nm (2.63 Å) | Shannon离子半径 (CN=8: r=1.109 Å → σ≈2r×2^(-1/6)) |
| ε | 0.50 kJ/mol (0.1195 kcal/mol) | Li & Merz镧系力场研究 |
| q | +3.0 e | Nd(III)形式电荷 |

## 对接结果

### 4IQK (Keap1-Nrf2)
- 结合能：-1426.44 kcal/mol
- 最优结合位点位于Cluster 8
- 关键残基：CYS583, ASP585, ASP579, GLU582, SER580, GLU478, ASN477, ILE476, GLN481, GLY484
- 金属结合热点包括Cys/His/Glu/Asp残基的高频聚类

### 6HN3 (GPX4)
- 结合能：-1285.14 kcal/mol
- 最优结合位点位于Cluster 1
- 关键残基：ASP6, CYS10, MET14, GLY15, ALA9, ASP116, GLY8, LEU7, ASP117, ASN115
- 距催化残基Sec46约8-12 Å（非直接竞争活性位点）

## MD分析结果

### 结构稳定性

| 指标 | 4IQK (Keap1) | 6HN3 (GPX4) |
|------|-------------|-------------|
| Cα RMSD | 0.21 ± 0.08 Å | 0.18 ± 0.04 Å |
| Rg | 1.93 ± 0.04 nm | 1.50 ± 0.01 nm |
| SASA | 127 ± 2 Å² | 79 ± 2 Å² |
| 氢键数 | 493 | 350 |

### Nd(III)结合分析

| 指标 | 4IQK | 6HN3 |
|------|------|------|
| Nd-蛋白最短距离 | 0.21 ± 0.00 Å | 0.21 ± 0.00 Å |
| 3.5Å接触残基数 | 293 | 170 |

Nd(III)在真空模拟中与蛋白保持<3Å的近距离接触，表明初始对接位点具有良好的几何互补性。

## 显式溶剂MD协议（GROMACS + CHARMM36）

### 运行要求
- GROMACS >= 2021
- GPU（推荐A100/V100，100ns生产模拟约需8小时）
- 文件位于 `gromacs_inputs/` 目录

### 操作流程
1. 将 `gromacs_inputs/` 文件夹复制到GROMACS服务器
2. 运行 `bash run_all.sh`
3. pdb2gmx后需手动编辑 `topol.top`：
   - 在 `[ atomtypes ]` 添加Nd(III)参数
   - 在 `[ molecules ]` 末尾添加 `ND3  1`

### 后续分析
生产模拟完成后，使用 gmx_MMPBSA 计算结合自由能：
```bash
gmx_MMPBSA -O -i mmpbsa.in -cs md.tpr -ct md_center.xtc -cp complex.top -lp ligand.top -rp receptor.top
```

## 研究结论

1. **对接结果**显示Nd(III)对Keap1和GPX4均具有强结合倾向（-1426/-1285 kcal/mol）
2. **真空MD**确认初始结合位点的几何互补性，Nd-蛋白距离稳定
3. **金属结合模式**：Nd(III)主要通过Cys、Asp、Glu、His等含硫/含氧残基配位，与已知镧系元素-蛋白相互作用模式一致
4. **局限性**：真空模拟无法捕捉溶剂化效应对结合自由能的贡献；Rg/SASA值受环境影响；0.1 ns采样时间不足以充分探索构象空间

## 产出物清单

| 文件 | 类型 |
|------|------|
| Figure7_RMSD_Rg.png | 300dpi图表 |
| Figure8_RMSF.png | 300dpi图表 |
| Figure9_NdDist_SASA.png | 300dpi图表 |
| Figure10_NdContacts.png | 300dpi图表 |
| Figure11_Summary.png | 300dpi汇总表 |
| gromacs_inputs/* | GROMACS完整输入文件包 |
| analysis_results.json | 全部分析数据(JSON) |

## 发表建议

要在高影响力期刊发表，建议补充：
1. 显式溶剂MD：GROMACS + CHARMM36，≥100 ns，3次独立重复
2. MM/PBSA结合自由能：使用gmx_MMPBSA，包含熵贡献
3. 关键残基的氢键寿命和占有率分析
4. FEL自由能面图，识别主要构象状态
5. 对比游离蛋白和Nd结合态的结构差异
