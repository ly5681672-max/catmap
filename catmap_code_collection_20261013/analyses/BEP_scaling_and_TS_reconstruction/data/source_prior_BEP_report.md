# O/H → TS1/TS2 → 环氧化 TOF：分析报告

## 结论先行

- 13 个有效位点均通过 `surface_name` 对齐；Ti–Zr 继续排除。输入温度为 333.15 K，TOF取自既有严格基准计算。本轮没有重跑稳态动力学、没有改能量/机理/条件。
- 旧表中的 O、Ha、Hb 只作为**历史参考描述符**使用，未做 slab 偏移修正。不同历史/严格模型参考仍未完成物理一致性验证，因此这里的映射和活性图均为探索性结果。
- 对 TS1 的直接 O/H 线性映射，最佳样本内组合是 O+Ha，但 `R²=0.138`、RMSE `0.632 eV`；LOOCV RMSE `0.780 eV`。这不是有力的 TS1 重建关系。
- 对 TS2 的直接映射，最佳组合是 Ha+Hb，`R²=0.695`、RMSE `0.429 eV`，但 LOOCV RMSE `0.754 eV`、LOOCV R² `0.057`，泛化稳定性有限。
- 实际 NEB 反应自由能→活化能的线性 BEP 中，NEB1 拟合较弱（R² `0.205`、RMSE `0.418 eV`）；NEB2 RMSE `0.117 eV`、R² `0.182`，低 RMSE 部分来自其活化能变化范围较窄，LOOCV R²仍为负值（`−0.148`）。不要只据 R² 或 RMSE 将其解释为经验证的普适 BEP。
- 将 O/H 映射到 TS1/TS2 后，经验二次表面和既有 CatMAP VectorMap 的最高值均落在 O/H 绘图区边界；**未发现被当前范围两侧数据支持的内部火山峰**。边界最大值不是火山峰。

## 数据与方法

分析数据来自归档的 `catmap_official_run/energies.txt`、`analyses/empirical_O_H_descriptor_analysis/data/matched_data.csv`、`reference_slab_comparison.csv`，以及 TS1–TS2 经验二次模型、已保存的 CatMAP 50×50 收敛网格及其官方 VectorMap 400×400 插值表。输入快照存于本目录 `data/`，复现脚本为 `scripts/analyze_oh_to_ts_to_tof.py`。归档中 `ts1_ts2_descriptor_grid.mkm` 是早期 5×5 pilot 配置模板；实际 50×50 运行由 `run_official_grid.py` 的网格参数生成，日志、网格摘要和运行脚本另存于 `catmap_official_run/O_H_to_TS1_TS2_to_TOF/`，不把 pilot 配置误称为 50×50 实际配置。

直接映射对 TS1_s、TS2_s 分别使用 O、Ha、Hb 的三种单描述符，及 O+Ha、O+Hb、Ha+Hb 三种双描述符；O+Ha+Hb 仅作三变量诊断。使用中心化/标准化最小二乘；样本内指标以 13 个催化剂为单位，LOOCV 每次只用训练折估计标准化量和回归系数。条件点不被当作独立催化剂重复计数。

串联路径定义：

1. NEB1/path2：`ΔG₂ = G(EpoxideHHO*) − G(Rpair*)`，`ΔG‡₂ = G(TS1*) − G(Rpair*)`；拟合 `ΔG‡₂ = αΔG₂ + β`。
2. NEB2/path4：`ΔG₄ = G(OHHa*) − G(HHO*)`，`ΔG‡₄ = G(TS2*) − G(HHO*)`；拟合 `ΔG‡₄ = αΔG₄ + β`。

链式预测先用 O/H 估计初态形成能与反应自由能，再用对应路径 BEP 估计 TS；LOOCV 中被留出的催化剂不使用其真实初态、反应自由能或 TS。负的原始 `TS−IS` 保留，不裁剪为正势垒：NEB1 中 Ti–W、Ti–Ta 共 2 个负值；NEB2 中 Ti–Mo、Ti–V、Ti–Cr 共 3 个负值。其物理解释需回到原始 NEB/结构核查，不能当作新势垒。

活性映射有两类，均是将 O/H 经线性模型映射至 TS1/TS2 后再取值：

- `empirical_TS_quadratic`：复用已存 TS1–TS2 经验二次响应面；不是 CatMAP 网格计算。
- `saved_CatMAP_VectorMap_interpolation`：查值于已有 CatMAP 求解网格的 VectorMap 插值；没有新增 CatMAP 求解。仅在既有 TS 描述符网格范围内评估，范围外保留空值，不边界夹取。

## 定量结果

### O/H 对 TS 态的直接预测

| 目标态 | 样本内最佳组合 | R² | RMSE (eV) | 最大绝对误差 (eV) | LOOCV RMSE (eV) | LOOCV R² |
|---|---:|---:|---:|---:|---:|---:|
| TS1_s | O+Ha | 0.138 | 0.632 | 1.388 | 0.780 | −0.313 |
| TS2_s | Ha+Hb | 0.695 | 0.429 | 1.048 | 0.754 | 0.057 |

这些数字只描述 13 个样本上的线性重建，不证明 O/H 是足以决定过渡态的物理坐标。特别是 Ha/Hb 可能彼此相关，不能把两者拟合系数直接解读为独立因果贡献。

### 实测反应自由能与过渡态的 BEP 诊断

| 路径 | α | β (eV) | R² | RMSE (eV) | LOOCV RMSE (eV) | LOOCV R² | 原始 TS−IS 负值 |
|---|---:|---:|---:|---:|---:|---:|---|
| NEB1/path2 | 0.367 | 0.907 | 0.205 | 0.418 | 0.713 | −1.318 | Ti–W、Ti–Ta |
| NEB2/path4 | 0.079 | 0.200 | 0.182 | 0.117 | 0.139 | −0.148 | Ti–Mo、Ti–V、Ti–Cr |

低样本量下这些线性关系不宜被提升为普适物理规律。NEB2 的拟合 RMSE 较小但 R²低，说明绝对误差较小与解释方差较弱可以同时成立。

### O/H→初态/反应能→BEP 的链式重建

| 路径 | 样本内表现较好的描述符 | R² | RMSE (eV) | LOOCV RMSE (eV) | LOOCV R² |
|---|---|---:|---:|---:|---:|
| NEB1/path2 | Ha（最佳单变量） | 0.128 | 0.636 | 0.706 | −0.075 |
| NEB2/path4 | Ha+Hb（最佳样本内组合） | 0.692 | 0.431 | 0.833 | −0.149 |

NEB2 的 O+Hb 链式模型样本内略逊（R² `0.689`、RMSE `0.433 eV`），但 LOOCV RMSE `0.546 eV`、LOOCV R² `0.506`；相比之下 Ha+Hb 的样本内排序更好而 LOOCV 较差。这一差异显示当前小样本的组合排名对评价目标敏感，不能只按样本内分数声称唯一最佳机制描述符。

### O/H 坐标下活性图的映射峰值

坐标边界均按 13 个历史 O/H 数值范围各外扩 5%。所有下列峰均为**边界最大值**，不是内部局部峰/火山峰。

| O/H 组合 | 活性背景来源 | 有效网格点 | 最大值坐标 (eV) | 最大 log10(TOF) | 最近真实位点* |
|---|---|---:|---|---:|---|
| O+Ha | TS 经验二次响应面 | 160,000/160,000 | O=2.244506, Ha=−3.388277 | −1.453217 | Ti–Co |
| O+Ha | 已存 CatMAP VectorMap 插值 | 160,000/160,000 | O=2.244506, Ha=−4.033044 | −1.448830 | Ti–Ni |
| O+Hb | TS 经验二次响应面 | 160,000/160,000 | O=1.456022, Hb=−4.722527 | 0.601332 | Ti（单 Ti） |
| O+Hb | 已存 CatMAP VectorMap 插值 | 155,508/160,000 | O=2.244506, Hb=−4.722527 | 4.823461 | Ti（单 Ti） |

\*最近点按各轴绘图范围归一化后的欧氏距离排序；仅表示几何邻近，不代表该催化剂达到预测活性，也不代表该区域经过新 CatMAP 求解验证。

两类图的峰位和峰高差异很大。例如 O+Hb 的既有 CatMAP 插值最大 log10(TOF)=4.823，而复用经验二次面仅为 0.601；这是两种响应模型在映射坐标上的差异，不能混写为同一个结果。O/H 到 TS 的直接重建误差偏大，加上历史 slab 参考未统一，使这些表面尤其不适合作为经验证的连续描述符预测。

## 参考能与科学限制

`reference_slab_comparison.csv`显示历史 O/H 形成能与严格 CatMAP 数据的参考可能有显著位点差别（表中若干位点差异达到数十分之一至约 1 eV 量级）。本轮未对 O/Ha/Hb 做平移或重算，因而不可把 O/H 旧参考与当前网络态参考视为已证明同源。对原始 Excel 的源单元格/公式和结构等价性，本轮也没有重新审计；当前证据只足以支持“按 surface_name 对齐的历史描述符探索”。

TS1/TS2 的能量来自当前 `energies.txt`；该输入及先前模型对其热修正/参考条件仍带有项目中记录的条件性限制。本次保持原表值，不再派生新能量或宣称共同参考已被独立验证。两条 NEB 的能量零点能否直接比较，也不由 BEP 相关性本身证明。

## 输出索引

- `data/linear_descriptor_model_metrics.csv`：O/Ha/Hb → 七网络态的单/双描述符回归评价。
- `data/linear_descriptor_model_predictions.csv`：逐催化剂样本内及 LOOCV 预测/残差。
- `data/BEP_model_summary.csv` 与 `data/BEP_catalyst_values_and_predictions.csv`：两路径 BEP 系数、误差、逐位点数据，保留负差值。
- `data/O_H_to_BEP_chain_model_metrics.csv`、`data/O_H_to_BEP_chain_TS_predictions.csv`：严格留一的链式重建评价及逐点结果。
- `data/mapped_activity_peaks.json`：四张映射图的范围、有效点数与峰值元数据。
- `data/network_state_energy_long.csv`、`data/matched_analysis_table.csv`：统一主表与网络态长表。
- `figures/scaling/`：直接态预测、路径 BEP、链式与实际 TS 比较图。
- `figures/tof_maps/`：O+Ha/O+Hb 的经验二次面和既有 CatMAP VectorMap 映射图，均含 13 个标签和边界最大值标记。
- `origin_exports/maps/`：400×400 XYZ、矩阵数据、13 位点坐标/TOF与峰值表；经验面 XYZ 全矩形无空值，CatMAP 插值在已有 TS 网格范围外留空。
- `scripts/analyze_oh_to_ts_to_tof.py`：可复现脚本；运行时只读本包快照并重新生成分析/绘图，不调用 `ReactionModel.run()`。

## 根因判断与建议

主要瓶颈不是绘图，而是历史 O/H 描述符并未显示对 NEB1 TS1 的足够预测信息；NEB2 的直接 TS2 样本内关系较好，但 LOOCV 明显下降。将 O/H→反应态→BEP 串联会累积初态、反应自由能和 BEP 三处误差，当前链式预测不能替代 TS 数据。当前最小下一步是先解决历史 O/H 与严格模型 slab/原子参考之间的可追溯关系，并审查负 TS−IS 的原始结构/能量定义；在此之前不建议把这些连续面作为论文定量预测或用作新催化剂筛选。此报告不启动新的 CatMAP 网格求解。

## 方法来源说明

本分析遵循可追溯数据处理、留一验证和可视化不夸大外推范围的流程。研究代理流程方法参考：Kassis et al., “Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents,” arXiv:2609.00065 (2026), https://doi.org/10.48550/arXiv.2609.00065。该文仅为分析流程参考，不是本项目机理或数据的科学证据。
