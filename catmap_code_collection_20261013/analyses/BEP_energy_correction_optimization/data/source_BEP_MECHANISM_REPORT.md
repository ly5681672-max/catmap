# BEP关系、共吸附修正与TS重建：结果报告

## 结论摘要

| 路径 | BEP关系（有符号） | R² | RMSE (eV) | LOOCV RMSE (eV) | LOOCV R² |
|---|---|---:|---:|---:|---:|
| NEB1 / Rpair→TS1→EpoxideHHO | `ΔE‡signed = 0.36719662 ΔG_rxn + 0.90724224` | 0.2054 | 0.4177 | 0.7134 | -1.3180 |
| NEB2 / HHO→TS2→OHHa | `ΔE‡signed = 0.07921318 ΔG_rxn + 0.19951166` | 0.1821 | 0.1172 | 0.1388 | -0.1483 |

`ΔG_rxn = G_FS − G_IS`，`ΔE‡signed = G_TS − G_IS`。这两个关系是13个既有位点上的样本内拟合；LOOCV表现较弱，尚不能作为可靠外推预测器。NEB1负的原始TS−IS值出现在：Ti-W (-0.566651 eV), Ti-Ta (-1.174636 eV)。NEB2负值出现在：Ti-Mo (-0.012820 eV), Ti-V (-0.024677 eV), Ti-Cr (-0.193155 eV)。这些是输入数据的有符号端点差，不被替换成零，也不被重新命名为正活化能。

## 单额外共变量的诊断

在基准BEP反应能之外逐次加入一个候选能量特征。当前样本内RMSE最低的候选：NEB1为`OH_star_legacy`（N=13，R²=0.4934，RMSE=0.3335 eV，LOOCV RMSE=0.7821 eV）；NEB2为`OH_star_legacy`（N=13，R²=0.4765，RMSE=0.0937 eV，LOOCV RMSE=0.1143 eV）。该候选为历史单吸附态OH*能量，当前未与严格模型slab统一；故仅是探索性相关，不据此宣布物理修正已成立。对所有特征的定义与参考状态限制见`origin_exports/BEP_extended_model_ranking.csv`。NEB1首步的`ΔG1=G(Rpair*)−G(H2O2_g)−G(C6H12_g)`是整体共吸附形成热力学差，不是纯粹的共吸附相互作用能。

## Ti–Fe敏感性与共吸附解释

去掉Ti–Fe后，NEB1基准R²=0.2063、RMSE=0.4345 eV；NEB2基准R²=0.1667、RMSE=0.1218 eV。该删除敏感性检查仅作诊断，正式主结果仍包括全部13个位点。严格Rpair态相对于气相反应物的ΔG1已导出；其解释包含各自吸附与共吸附整体效应，不能独立拆成“相互作用能”。

## TS重建与既有CatMAP表面映射

`BEP_TS1_TS2_reconstruction.csv`分别保留使用全样本拟合和LOOCV barrier预测、再加回真实对应IS能量的TS重建结果。前者是in-sample诊断，后者是更保守的留一预测。另复用既有O/H描述符链的LOOCV预测（NEB1 Ha、NEB2 Hb），按其TS1/TS2坐标查询既有400×400 VectorMap插值表：13/13个位点坐标落在保存网格范围内。导出的TOF仅为最近插值网格点值，不是新CatMAP求解，也不是重建TS的真实TOF；网格外不赋值。底图的极值或网格质量仍以原有50×50直接求解审计为准。

## 科学判断与限制

1. 两条传统BEP线性关系在样本内解释度偏低，LOOCV不支持稳健跨位点预测。不能把全局线性BEP当作已验证普适定律。
2. NEB1样本内拟合误差明显较大；NEB2绝对误差较小但LOOCV仍有限。负的TS−IS端点差保留，可能反映局部端点能量/构型排序，不应被端点差直接解释成通常意义的正活化势垒。
3. OH*额外项在两条反应的样本内拟合均有改善。NEB1的LOOCV RMSE由0.7134增至0.7821 eV，未改善；NEB2则由0.1388降至0.1143 eV，呈现有限的留一改善。但OH*历史slab参考闭合不完整，样本仅13个位点，因此尚不能升级为正式共吸附修正或外部有效的BEP模型。
4. 旧O/H单物种数据与严格六步态能量存在位点相关参考偏移；尤其Ti–Fe历史clean-slab差异约1.064927 eV（历史审计记录），当前缺少逐原子结构、cell和同一热修正口径的证据来解释该差值。因此旧吸附能仅作探索性变量，不能与严格Rpair/TS绝对能量直接混成统一参考。
5. 首步共吸附状态的真实结构与独立吸附态简单叠加并不等价。模型中的Rpair自由能可以作为整体状态量，但当前资料不足以单独分离吸附能和吸附物间相互作用能。
6. 官方CatMAP B.E.P. scaler系数已在隔离模式下复现；它验证的是实现/拟合系数，不验证拟合质量或物理普适性。CatMAP文档所述BEP是TS能相对初态对反应能的标度关系，具体解释仍依赖一致的参考能和物种定义。

## 复现与产物

- `scripts/verify_catmap_bep.py`：仅调用官方Scaler系数矩阵，无`ReactionModel.run()`。
- `scripts/analyze_bep_reconstruction.py`：重建13位点BEP、扩展项、留一评估、TS重建和既有VectorMap近邻读取。
- `scripts/render_bep_figures.py`：科研图及Origin表输出。
- `data/official_BEP_scaler_metadata.json`：注明实际调用状态。
- `reports/SHA256SUMS.csv`：本归档包文件哈希。

未做：新的DFT/NEB、CatMAP稳态/TOF、反应条件搜索、CatMAP网格计算。
