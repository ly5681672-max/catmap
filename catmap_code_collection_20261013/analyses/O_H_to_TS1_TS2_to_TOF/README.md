# O/H 到 TS1/TS2 重建分析包

本包复用 13 个严格基准收敛位点、历史 O/Ha/Hb 描述符、既有 NEB 能量及既有 CatMAP TS1–TS2 活性网格；不包含新的 DFT、稳态求解或条件搜索。

- `data/`：只读输入快照、surface_name 匹配表、模型指标、逐位点预测、网格映射峰值。
- `figures/`：直接标度、BEP/链式诊断及 O+Ha/O+Hb 活性映射图。
- `origin_exports/maps/`：400×400 XYZ、矩阵、真实点及峰值 CSV。
- `scripts/`：分析与 SHA256 清单生成脚本。
- `reports/`：中文结果与科学限制。
- `catmap_official_run/`：既有 CatMAP 官方 50×50 网格的脚本、日志、输入和摘要快照；不是本轮新运行。
- `SHA256SUMS.csv`：本包文件及原始 H 盘来源文件的 SHA256。

核心解释限制：O/Ha/Hb 仍为历史参考，未做 slab 平移；任何二维颜色面均为探索性映射，不能当作已验证的 O/H 连续 CatMAP求解结果或内部火山峰。
