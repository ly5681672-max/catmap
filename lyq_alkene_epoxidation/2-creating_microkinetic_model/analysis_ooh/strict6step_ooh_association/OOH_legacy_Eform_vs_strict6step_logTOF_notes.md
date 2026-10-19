# 历史 OOH* 形成能—正式六步 TOF 关联（探索性）

- N = 13；仅按 `surface_name` 一对一匹配，13/13 位点一致。
- 描述符：`OOH_star_legacy`，来源为 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization\data\source_carbon_O_H_TS1.csv`，参考字段均为 `DFT_Gibbs_single`。
- 目标：正式六步 CatMAP 结果中 `C6H12O_g` 的 turnover frequency；输入表 13/13 状态为 CONVERGED、333.15 K、正向。TOF 与输入表内 log10 值最大差异小于 5e-16。
- Pearson r = 0.548421, p = 0.052301; Spearman rho = 0.543956, p = 0.054647.
- 样本内直线描述：R2 = 0.300765, RMSE = 6.925243 log10 units, MAE = 4.985862 log10 units；这只是简要描述统计，不是机理拟合或因果验证。
- 绝对线性残差较大的样本（仅作诊断线索）：Ti-Hf, Ti-Re, Ti-Ni, Ti。

## 解释边界

OOH* 是历史独立吸附态能量描述符；严格六步网络不含 OOH* 物种，因此该图不提供六步模型中的 OOH* 覆盖度。两边能量/TOF来自不同模型与参考数据路径，尚未证明slab与热修正完全统一。图中相关性只能用于探索，不能解释为八步 OOH* 稳态覆盖度导致 TOF 的因果证据。

## 文件

- 配对数据：`H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\strict6step_ooh_association\data\OOH_legacy_Eform_vs_strict6step_logTOF_13sites.csv`
- 图：`H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\strict6step_ooh_association\figures\OOH_legacy_Eform_vs_strict6step_logTOF_13sites.png`、`H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\strict6step_ooh_association\figures\OOH_legacy_Eform_vs_strict6step_logTOF_13sites.pdf`
- 复现脚本：`H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\strict6step_ooh_association\scripts\plot_ooh_legacy_vs_strict6step_tof.py`
