# TS1* 的 C/O/H 双能量描述符探索（333.15 K 项目条件）

## 结论摘要

- **独立 C*：未找到。** 两份 Excel 的所有工作表、旧 `DFT_formation_catmap.txt`、转换脚本及归档中可检索的计算文件名都没有提供 13 个位点的原子 C 吸附态及其可追溯总能。因此 C、C+O、C+Ha、C+Hb 的数值拟合均未执行；没有以 C6H12* 或 μC 替代 C*。
- 在可计算的三组 O/H 双描述符中，`O+Ha` 的样本内 R² 最高且 RMSE 最低：**R²=0.1382、RMSE=0.6321 eV、MAE=0.5281 eV、LOOCV RMSE=0.7803 eV（N=13）**。它相对 Ha 单描述符（R²=0.1355、RMSE=0.6330 eV）几乎没有改善；不构成强 TS1 描述符关系。
- 因 C* 数据不存在，无法回答“加入 C* 是否改善 O+Ha”，也不能在六种双描述符中判定真正全局最佳。当前**最佳可计算候选**是 O+Ha，而不是六者中的已证最优。
- 未生成 C 相关响应面或 C 描述符 Origin 网格；这会把未知 C* 值伪装成数据。已输出全部可用 O/H 线性模型的比较图、催化剂命名的实际-预测图和 Origin 点/系数/残差表。

## 数据和定义

- 目标 `TS1_formation_energy_eV` 取自归档 `catmap_official_run/energies.txt` 中 13 个已纳入位点的 `TS1` 行；`tizr` 明确排除。历史 `O`、`Ha`、`Hb` 与 TOF 来自 `empirical_O_H_descriptor_analysis/data/matched_data.csv`，按 `surface_name` 一对一连接，不按行序连接。匹配数据标记 13/13 `CONVERGED`、333.15 K。
- 回归形式严格为普通线性最小二乘：单变量 `TS1=aX+c`；双变量 `TS1=aX+bY+c`。没有高阶项、变量筛选或异常点剔除。LOOCV 每折只用其余 12 个位点重拟合。
- O/Ha/Hb 是旧的单物种历史形成能，不是当前六步网络内的独立反应态，也不等于 `EpoxideHHO`/`HHO`。历史形成能表公式可重建并与旧表相符（39/39），但历史单物种的温度记录 39/39 缺失。严格模型 TS1 与历史 O/H 的 Gibbs 温度/结构参考不能据此视为完全同源。
- 参考 slab 审计显示旧/当前 slab 总能差并非统一零偏移：最大示例 Ti-Fe `+1.064927 eV`、Ti-Co `+0.877702 eV`、Ti-Ti 约 `+0.000005 eV`、Ti-Ta `+0.128845 eV`、Ti-Ni `−0.114951 eV`。这些是审计中记录的 `strict−legacy` slab 数值差，不能直接当作描述符校正量；本分析不对 O/H 或 TS1 作任何能量平移。
- 转换脚本将 `μC=−9.28 eV` 写成“原始 DFT 原子参考能”，CatMAP energies reference 也沿用该常量；但未找到相应孤立 C/参考碳相计算输出或自旋、泛函等设置证据。该常数不是 C* 能量，当前不足以为 C* 形成能提供完整参考证明。

## 可用回归结果

系数均以表中 eV 描述符直接拟合，目标为 eV。

| 模型 | N | 斜率/系数 | 截距 eV | 样本内 R² | RMSE eV | MAE eV | 最大绝对误差 eV | LOOCV RMSE eV |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| O | 13 | O: −0.017964 | −33.268616 | 0.0020 | 0.6802 | 0.5731 | 1.3324 | 0.7982 |
| Ha | 13 | Ha: +0.275110 | −32.805991 | 0.1355 | 0.6330 | 0.5260 | 1.3580 | 0.7150 |
| Hb | 13 | Hb: +0.004175 | −33.253277 | 0.0001 | 0.6808 | 0.5734 | 1.3131 | 0.8438 |
| O+Ha | 13 | O: −0.020621；Ha: +0.275806 | −32.814726 | **0.1382** | **0.6321** | 0.5281 | 1.3880 | 0.7803 |
| O+Hb | 13 | O: −0.026128；Hb: −0.015890 | −33.298142 | 0.0027 | 0.6799 | 0.5750 | 1.3185 | 0.8999 |
| Ha+Hb | 13 | Ha: +0.279285；Hb: −0.020317 | −32.831849 | 0.1372 | 0.6324 | **0.5224** | 1.3260 | 0.9381 |

`C`、`C+O`、`C+Ha`、`C+Hb` 的 N=0/系数/R²/RMSE 均为“不可估计”，不是零；availability 表记录为 `NOT_FIT_CSTAR_MISSING`。三组已有 O/H 对照模型使用完全相同的 13 个位点。

### 当前可回答的比较

- 样本内 R²/RMSE 最佳可计算双描述符：O+Ha，但与 Ha 单变量几乎等同（R²仅增加 0.0027，RMSE 降低约 0.0010 eV）。这说明在当前 13 点、历史 O/H 参考数据上，没有证据表明 O 对 Ha 提供了明显额外的 TS1 重建能力。
- Ha+Hb 的 MAE 略低，但 R²和 RMSE略差，LOOCV RMSE更高；按“样本内 R²及 RMSE”为主排序，仍列 O+Ha 第一。该差异很小，不应过度解读。
- 各模型最大误差约 1.31–1.39 eV，LOOCV RMSE约 0.71–0.94 eV。样本内关联整体弱；并且描述符形成能温度/参考不完整，当前结果只用于筛选下一步数据价值，不构成可外推的物理标度关系。

## C* 缺失时最少要补什么

对每个目标位点，需要一个明确单原子 C* 吸附构型的 `G(slab+C*)`，以及严格同晶胞、同计算参数/热修正的 `G(clean slab)`。若已有 clean slab 被证明同源可复用，则新增 13 个 slab+C* 计算即可；若不能证明，则 13 个 clean slab 也需提供/重算。还要给出所用 `μC` 的可追溯参考（来源相、计算方法、自旋等）及统一温度修正，才能计算 `G_form(C*)=G(slab+C*)−G(clean slab)−μC`。详细字段见空白模板。

## 输出索引

- C* 源扫描：`data/Cstar_source_inventory.csv`
- 13 位点 C* 原始数据模板：`data/Cstar_13site_DFT_data_template.csv`
- 匹配后的 TS1/O/H/TOF 表：`data/matched_TS1_O_H_TOF_data.csv`
- 全模型指标和逐点残差：`data/descriptor_fit_metrics.csv`、`data/descriptor_fit_predictions.csv`
- 可用模型实际-预测图：`figures/O_Ha_TS1_actual_vs_predicted.png/.pdf`、`figures/O_Hb_TS1_actual_vs_predicted.png/.pdf`、`figures/Ha_Hb_TS1_actual_vs_predicted.png/.pdf`
- 性能比较图：`figures/descriptor_pair_fit_performance.png/.pdf`
- Origin 文件：`origin_exports/catalyst_points.csv`、`origin_exports/fit_coefficients_and_metrics.csv`、`origin_exports/predictions_residuals.csv`
- 所有输入/输出 SHA256：`data/SHA256SUMS.csv`
- C* 源证据与参考限制：`reports/CSTAR_SOURCE_AUDIT.md`

未生成 C 相关 200×200/400×400 TS1 响应面，因为没有任何真实 C* 坐标；也没有运行 CatMAP、重算 TOF 或开展新 DFT。

## 方法学备注

小样本模型报告同时保留样本内拟合和逐点留一预测，避免把拟合优度等同于外推能力。该分析没有进行六种模型选择后的显著性推断或多重搜索校正。方法学参考：Kassis et al., *When does model selection outperform model averaging?* (2026), [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065)。此引用只用于小样本验证设计背景，不是本项目化学结论的数据来源。
