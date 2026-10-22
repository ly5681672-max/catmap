# 科研模型约束协议

- `SCIENTIFIC_MODEL_BASELINE.json`：冻结的正式八步输入和反应基准。
- `scripts/validate_scientific_model.py`：运行前检查指定 `.mkm` 路径、冻结哈希及 CatMAP 0.3.1 解析的反应、物种、`n_sites`、温度和储库字段；在临时目录中解析，不修改科研输入。
- `tests/test_validate_scientific_model.py`：验证基线通过、模型/能量/储库变更及非正式模型路径会被拦截。

每次计划运行时，将实际 `.mkm` 路径传给预检器。正式模型变更只能在研究者明确批准后登记为独立版本。不要为使预检通过而刷新当前冻结基准。
