# 本轮分析包索引

范围：复核现有CatMAP 0.3.1严格基线和既有TS1–TS2官方网格缓存；温度333.15 K；13个已收敛位点；Ti–Zr排除。

## 核心报告

- `EXISTING_WORK_REUSE_AUDIT.md`：既有结果复用/纠错分类。
- `TIFE_DATA_AND_MECHANISM_AUDIT.md`：Ti–Fe原始G、参考slab、模型换算和已有稳态/DRC核对。
- `MECHANISM_CORRELATION_REPORT.md`：ΔG3、覆盖度、DRC、TOF与缓存峰两侧的合并结论及证据边界。

## 关键数据

- `../data/combined_dG3_coverage_DRC_TOF.csv`：13位点严格基准合并表；按surface_name一对一匹配。空位来自高精度`coverage_site_balance.csv`。
- `../data/TIFE_13site_vs_12site_comparison.csv`：完整13点并含12点敏感性标记；不是剔除Ti–Fe后的替代基线。
- `../data/TIFE_12site_sensitivity_only.csv`：仅供敏感性比较。
- `../data/tife_workbook_cell_audit.csv`：两版工作簿哈希、Ti–Fe公式/缓存G、温度和路径单元格。
- `../data/TS1_TS2_peak_mechanism_data.csv`：缓存峰及轴向邻点的CatMAP官方Scaler网络态能量、有效ΔG2/3/4与显式TS步骤正向势垒。
- `../data/official_scaler_anchor_reproduction.csv`：13个位点、六步官方能量API与严格基线的交叉核验。
- `../data/source_manifest_sha256.csv`：原始输入/缓存的SHA256。

## 图和Origin表

`../figures/`含PNG/PDF：13位点ΔG3-覆盖度-TOF联合关系、覆盖度/DRC、TS1/TS2峰两侧覆盖度与TOF、峰两侧有效反应能、Ti–Fe代表位点比较。

`../origin_exports/`含合并点表、峰值及两侧覆盖度/TOF、官方Scaler能量与反应能、2500节点既有覆盖度缓存。缓存节点是此前CatMAP计算，不是本轮新求解。

## 复现和限制

- `../scripts/build_mechanism_package.py`只读已有CSV/Excel/energies.txt/官方pickle，调用CatMAP 0.3.1官方Scaler及ReactionModel反应能函数；不调用`run()`、求解器、Mapper或条件搜索。
- `run_validation.log`记录本轮关键通过门控。
- Ti–Fe旧slab与严格slab相差1.064927 eV，但缺少结构输出证明二者等价；旧比较表将状态10的F41误标成alternate slab。本包更正解释，未改旧文件或能量。
- NEB2相关Excel温度单元格C24/C28/C36为空，本轮没有填造温度或热修正。
