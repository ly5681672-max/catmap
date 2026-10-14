# CatMAP 官方网格证据快照

本轮没有执行新的 CatMAP 动力学求解。此目录仅收录既有 TS1/TS2 官方 CatMAP 0.3.1 网格的可追溯材料，供分析复现：

- `run_official_grid.py`：既有网格运行入口；
- `catmap_ts1_ts2_grid_50x50.log`：50×50 实际运行日志；
- `grid_summary.json`：网格节点、稳态残差和配置摘要；
- `energies.txt`：当时使用的能量表快照。

新分析脚本读取 400×400 的既有 CatMAP `VectorMap` 插值表（位于 `../data/`），并将 O/H 通过本轮拟合的 TS1/TS2 映射后查询。此过程不是新的 CatMAP `model.run()`；网格外点保留空值。旧 O/Ha/Hb 参考尚未完成与严格模型的物理 slab 参考一致性验证。
