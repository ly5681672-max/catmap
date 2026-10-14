# 已有成果复用与复核审计

本轮只读复用CatMAP 0.3.1严格基准结果；未运行`ReactionModel.run()`、稳态求解、官方描述符网格或储库条件搜索。

| 任务项 | 状态 | 证据/处理 |
|---|---|---|
| 13位点 ΔG3–TOF | `VERIFIED_REUSE` | 由严格基准形成能独立重算标准ΔG3，并与已有step-3有效ΔG按surface_name逐一合并；旧的activity_analysis_round1关系已存在，本轮只作一致性复核。 |
| EpoxideHHO*/HHO*覆盖度 | `VERIFIED_REUSE` | 复用高精度coverage_13sites.csv；空位另取coverage_site_balance.csv，避免Ti–Fe被旧宽表舍入成0。 |
| 13位点官方DRC | `VERIFIED_REUSE` | 复用rate_control_matrix.csv；没有将DRC解释成唯一决速步。 |
| TS1–TS2峰值及四个轴向邻点 | `VERIFIED_REUSE` | 复用已存5点表及50×50 CatMAP缓存的覆盖度/TOF；未重新求解网格。 |
| 网格点ΔG2/ΔG3/ΔG4和正向势垒 | `NOT_DONE -> 本轮补齐` | 通过CatMAP 0.3.1 `GeneralizedLinearScaler.get_free_energies`与`ReactionModel.get_rxn_energy`对既有5个缓存坐标作能量评估，并用13个真实锚点与基线输出交叉核验。 |
| Ti–Fe clean slab历史/严格参考比较 | `NEEDS_CORRECTION` | 旧reference_slab_comparison.csv将F41列作alternate slab；实际F41是状态10（C6H12O+H2O），历史slab应为catmap数据集!G51。已在新审计中更正，不改旧文件。 |
| Ti–Fe历史与严格slab结构等价性 | `UNCERTAIN` | Excel路径可读，但本归档没有相应POSCAR/CONTCAR/OUTCAR来证明结构/组成/计算设置一致；不作数值校正。 |

数据来源和各文件SHA256见`data/source_manifest_sha256.csv`。输出不替换旧图或历史审计文件。
