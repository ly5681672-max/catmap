# 含碳吸附物种与 O/H 特征对 TS1 / NEB1 的探索性标度分析

## 范围与来源

- 13个位点：Ti-Fe, Ti-Mn, Ti-Hf, Ti-Re, Ti-Mo, Ti-V, Ti-Co, Ti-Ti, Ti-W, Ti-Ta, Ti-Cr, Ti, Ti-Ni；Ti–Zr排除。
- 使用归档的原始形成能表、当前严格模型 `energies.txt`、已有13位点收敛TOF表和slab对照审计；所有连接均按 `surface_name` 一对一完成。
- 原始形成能表与归档副本SHA256一致：`9ea18aa87759d3df9b3c268eebb5631f55fc71d4737605ab3b2d991cc4014baf`；原始转换脚本SHA256：`8bc8fa564a39e15498284e0c64b40417dc1115a18101c997419f6506d7bd1ba9`。Excel源单元格和公式见 `data/source_workbook_crosswalk.csv`。
- 没有重新运行CatMAP、稳态、条件搜索或DFT；所有图为普通线性回归的经验响应，不是CatMAP描述符网格。
- 单物种 `C6H12*` 与 `C*` 不同。本项目既有数据是完整1-己烯吸附态 `C6H12*`，并非单原子碳。

## 主要结果

- O+Ha历史参考基准复核：R²=0.1382, RMSE=0.6321 eV (N=13);与此前0.1382基准相符。
- TS1形成能的数值最优双特征（含网络态/所有主候选）：`C6H12O_star_legacy + OH_star_legacy` (N=13, R²=0.2586, RMSE=0.5862 eV, MAE=0.4456 eV, LOOCV RMSE=0.7342 eV)。若组合含 `Rpair*`，因目标 `TS1` 与反应初态高度相关，必须视为网络能量重建诊断，不自动视作简化独立描述符。
- TS1形成能的最佳“烯烃吸附态 + 反应前O/H稳定态”组合：`C6H12_star_legacy + Ha_star_legacy` (N=13, R²=0.1451, RMSE=0.6295 eV, MAE=0.4990 eV, LOOCV RMSE=3.0107 eV)；相对O+Ha样本内R²变化 +0.0069。
- 仅作参考敏感性诊断的 slab 平移轨迹中，TS1双特征最佳为：`C6H12O_star_slabshift_conditional + Ha_star_slabshift_conditional` (N=13, R²=0.5240, RMSE=0.4697 eV, MAE=0.3351 eV, LOOCV RMSE=0.6188 eV)；其中前驱C6H12*与O/H组合为 `C6H12_star_slabshift_conditional + Ha_star_slabshift_conditional` (N=13, R²=0.4093, RMSE=0.5233 eV, MAE=0.3860 eV, LOOCV RMSE=1.1693 eV)。这类提升不能在参考与热修正闭合前当作已验证规律。
- NEB1原始有符号差值目标 `TS1−Rpair` 的最佳双描述符：`Ha_star_legacy + OH_star_legacy` (N=13, R²=0.6276, RMSE=0.2859 eV, MAE=0.2233 eV, LOOCV RMSE=0.3916 eV)。
- 同一条件slab平移敏感性轨迹的NEB1差值最佳双特征为：`Ha_star_slabshift_conditional + OH_star_slabshift_conditional` (N=13, R²=0.6606, RMSE=0.2730 eV, MAE=0.2142 eV, LOOCV RMSE=0.3816 eV)；仅为参考敏感性比较。
- 仅在不把 `Rpair` 直接用于构造差值目标的含碳前驱候选中，NEB1差值最佳组合：`C6H12_star_legacy + OH_star_legacy` (N=13, R²=0.5175, RMSE=0.3255 eV, MAE=0.2485 eV, LOOCV RMSE=6.5545 eV)。负差值原样保留，没有把能垒截断为零。
- NEB1反应能对有符号 `TS1−Rpair` 的单变量BEP式拟合：R²=0.2054、RMSE=0.4177 eV、LOOCV RMSE=0.7134 eV；这是同一网络能量内部的关系，不等于独立验证。
- 若将 `Rpair*` 纳入预测变量，数值最佳网络态耦合诊断为TS1：`OH_star_legacy + Rpair_strict` (N=13, R²=0.7507, RMSE=0.3399 eV, MAE=0.2892 eV, LOOCV RMSE=0.4886 eV)；NEB1差值：`OH_star_legacy + Rpair_strict` (N=13, R²=0.4737, RMSE=0.3399 eV, MAE=0.2892 eV, LOOCV RMSE=0.4886 eV)。这些不列为简化的独立描述符结论，因为目标本身与该网络态直接相关。
- 最佳单一历史吸附态特征TS1模型（排除Rpair网络态）：`Ha_star_legacy` (N=13, R²=0.1355, RMSE=0.6330 eV, MAE=0.5260 eV, LOOCV RMSE=0.7150 eV)。`Rpair*`单特征单列为网络耦合诊断：`Rpair_strict` (N=13, R²=0.5586, RMSE=0.4524 eV, MAE=0.3271 eV, LOOCV RMSE=0.5076 eV)。

### C6H12* + O/H 对TS1形成能的逐组合结果（N=13）

| 组合 | 历史参考 R² | 历史 RMSE (eV) | 历史 LOOCV RMSE (eV) | 条件slab平移 R² | 条件slab平移 RMSE (eV) | 条件slab平移 LOOCV RMSE (eV) |
|---|---:|---:|---:|---:|---:|---:|
| C6H12* + O* | 0.0086 | 0.6779 | 4.1903 | 0.0764 | 0.6543 | 1.7091 |
| C6H12* + Ha* | 0.1451 | 0.6295 | 3.0107 | 0.4093 | 0.5233 | 1.1693 |
| C6H12* + Hb* | 0.0030 | 0.6798 | 4.1791 | 0.1193 | 0.6390 | 1.7352 |
| C6H12* + OOH* | 0.0099 | 0.6775 | 3.9041 | 0.1095 | 0.6425 | 1.9865 |
| C6H12* + OH* | 0.0278 | 0.6713 | 3.6390 | 0.1865 | 0.6141 | 1.6710 |

条件平移列仅用于显示参考选择敏感性；不是已验证的统一物理参考。

## 物种可用性与定义

`DFT_formation_catmap.txt` 中存在13个共同位点的 `C6H12*`、`C6H12O*`、`O*`、`Ha*`、`Hb*`、`OOH*`、`OH*` 和候选 `OOH-C6H12*`、`Hb-O*`。C6H12和C6H12O气相形成能也存在。`C6H12*`吸附自由能按 `Eform(C6H12*)−Eform(C6H12_g)` 得到；所用历史气相值分别为 -26.031223 eV 和 -27.667771 eV。气相项为跨位点常数，因此含截距的OLS中它仅平移截距，不改变R²/残差。没有独立C*，也未将 `C6H12*` 冒充C*。

## 能量参考与解释限制

旧单吸附态和当前六步网络能量并非已完全证明为同一物理参考。旧/当前洁净slab差值在13个位点间并不统一，已知大偏移包括Ti-Fe、Ti-Co、Ti-Ti、Ti-Ta、Ti-Ni及Ti-Mo。报告同时保留：

1. `legacy_reference`：直接使用历史表原数值；
2. `slabshift_conditional`：仅按已审计的当前减历史slab能差平移旧吸附态形成能，假设对应吸附态总G及热修正在其他方面可比。

第二种按 `Eform_conditional = Eform_legacy - (Gslab_current - Gslab_legacy)` 计算，仅为条件性参考对齐诊断，不代表热修正/结构等价性已经验证。所有排名同时保存两条轨迹。形成能的绝对目标 `TS1` 来自当前严格模型输入，仍沿用其原有条件性来源标记。

`OOH-C6H12*` 和 `Hb-O*` 在旧转换脚本中被标为 `DFT_Gibbs_TS_candidate`，注释明确其并不等同于共吸附NEB TS自由能。因其与TS1的结构身份/独立性未核实，OOH-C6H12只作为单变量信息重叠诊断，未进入主组合排名；Hb-O不作为稳定态描述符。`Rpair*`是当前网络前驱态，主排名中单独标识，并避免把它作为NEB1差值候选组合的机制结论。

## 图与Origin数据

- `figures/TS1_top_descriptor_performance.png/.pdf`：TS1形成能最佳候选模型指标。
- `figures/TS1_actual_vs_predicted_carbon_pair_*.png/.pdf`：最佳含碳前驱态对照模型，带13个位点标签及y=x线。
- `figures/TS1_actual_vs_predicted_best_carbon_pair_*.png/.pdf`：全体非循环含碳组合中的高排名对照；产品吸附态模型单独标作事后描述。
- `figures/NEB1_barrier_actual_vs_predicted_carbon_pair_*.png/.pdf`：NEB1原始有符号能量差拟合图。
- `figures/*_200x200.png/.pdf`、`*_400x400.png/.pdf`：经验线性响应面，无凸包遮罩、全矩形显示、标出全部13个数据点和名称；不是CatMAP连续网格。
- `origin_exports/*_xyz.csv` 与 `*_matrix.csv`：完整规则网格；`*_catalyst_points.csv`：真实点坐标、目标值与TOF。
- `data/fit_metrics_all_candidates.csv`、系数、逐位点预测/残差CSV和两个目标的排名CSV保存全部可用模型；`data/descriptor_source_inventory.csv`记录物种定义与来源。

## 科学解释边界

这些是13个催化剂的小样本内线性重建结果。LOOCV用于稳定性提示，不作为绘图门槛。拟合相关不证明因果机制；含产品吸附态的模型属于事后描述，含 `Rpair` 的模型对反应能垒目标有代数耦合。最优关联应结合结构核查、统一热力学参考和独立位点数据再用于TS1–TS2–TOF简化模型。当前没有证据把任何经验线性响应面称为火山峰。

## 可复现性

脚本：`scripts/analyze_carbon_ts1_scaling.py`。本分析没有重算任何TOF或CatMAP稳态。原始输入SHA256列于 `data/analysis_summary.json`；输出SHA256列于 `reports/SHA256SUMS.csv`。
