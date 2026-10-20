# 显式 OOH 八步模型验证与覆盖度筛选

本目录仅服务于 `explicit_OOH_8step` 研究主线。它复核原始 Excel/NEB 能量来源、使用本机 CatMAP 0.3.1 解析原始八步模型，并对已有的 Ti–Co 零库存诊断根进行**展示层微量覆盖度筛选**。不运行 CatMAP 稳态、TOF、DRC 或批量计算；不修改原始输入和历史结果。

## 目录

- `reports/`：模型状态、NEB 映射、位点守恒与计算就绪门控。
- `data/eight_step_coverage_screening.csv`：Ti–Co provisional 根逐稳定物种筛选表。
- `data/eight_step_parse_screen_summary.json`：CatMAP 解析、矩阵不变量和筛选汇总。
- `data/source_audit/`：此前实际 Excel 审计的只读副本，用于让本目录自包含；原始审计文件不变。
- `scripts/screen_tico_zero_inventory_coverages.py`：只解析模型并读已保存根，不调用 `.run()`。
- `tests/test_trace_coverage_rules.py`：筛选规则单元测试。

运行示例（在项目根目录，使用 CatMAP 0.3.1 环境）：

```powershell
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\scripts\screen_tico_zero_inventory_coverages.py `
  --model-dir .\lyq_alkene_epoxidation\2-creating_microkinetic_model `
  --root-json .\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\zero_inventory_tico_diagnostic\data\tico_zero_inventory_root.json `
  --output .\analyses\explicit_OOH_8step_validation\data\eight_step_coverage_screening.csv
```

该脚本的筛选输出只说明给定已保存根上的数值覆盖度与展示近似；不证明 NEB 物理映射、洁净表面全局可达性或跨催化剂普适性。
