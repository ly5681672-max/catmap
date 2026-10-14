# 既有BEP工作与数据入口审计

## 本轮复用的数据

- 严格CatMAP形成能表快照：`data/source_strict_energies.txt`，按`surface_name`与13个`CONVERGED`位点匹配；Ti–Zr排除。
- 13位点TOF及历史O/Ha/Hb：`data/source_matched_O_H_TOF.csv`。
- 既有BEP逐点值、回归摘要和拟合链预测：`data/source_prior_BEP_points.csv`、`data/source_prior_BEP_summary.csv`、`data/source_prior_chain_predictions.csv`。
- 含碳/单吸附态候选及参考说明：`data/source_carbon_O_H_TS1.csv`、`data/source_descriptor_inventory.csv`。
- 已有CatMAP 50×50直接求解网格和400×400 VectorMap插值：`data/source_catmap_grid_50x50.csv`、`data/source_catmap_VectorMap_400x400.csv`，连同`data/source_grid_summary.json`。
- 官方Scaler隔离核验产物：`data/official_BEP_scaler_coefficients.csv`与`data/official_BEP_scaler_metadata.json`。

## 数据流核验

13个位点由`surface_name`键连接，不用行号拼接。BEP输入由严格表中反应物/产物/TS形成能重新计算，并与既有逐点表核对。NEB有符号TS−IS差保留负值。历史单吸附态数据只作为探索性扩展特征，参考闭合未验证，不能解释成已统一到严格模型slab。

## 范围与禁止操作

本包不调用`ReactionModel.run()`，不重跑稳态、不扫条件、不更改DFT、NEB、TOF、反应机理或历史网格。CatMAP连续面仅复用已存VectorMap插值；不是新求解的CatMAP网格。本轮代码与输出先在C盘隔离目录生成，完成检查后再作为新目录复制入H盘归档；原始数据不覆盖。
