# 显式 OOH 八步模型：失败位点与物理映射诊断

此目录为独立、可复核的诊断快照，不替换 `multisite_batch_8step` 或正式模型。

## 结构

- `inputs/`：原始八步模型、候选模型及失败位点输入的只读快照。
- `data/`：签名 NEB 端点分类、CatMAP 解析矩阵/守恒零空间、网络平衡、覆盖度汇总、失败诊断、来源口径和 SHA-256 清单。
- `scripts/`：由 CatMAP 0.3.1 解析源模型并生成诊断表的脚本；Ti–Zr 有界初猜恢复尝试脚本；包清单生成器。
- `tests/`：对指定七条异常 NEB 及两种网络库存维数的回归测试。
- `logs/`：有限恢复尝试日志；不收敛/中断记录不是结果通过证明。
- `runs/`：隔离的候选模型运行输入与失败位点副本。
- `reports/`：汇总解释和物理门控报告。

## 复现

在项目根目录执行：

```powershell
G:\anaconda\python.exe analyses\explicit_OOH_8step_validation\multisite_failure_and_full_model_diagnostics\scripts\build_failure_physical_audit.py
G:\anaconda\python.exe -m unittest discover -s analyses\explicit_OOH_8step_validation\multisite_failure_and_full_model_diagnostics\tests -v
G:\anaconda\python.exe analyses\explicit_OOH_8step_validation\multisite_failure_and_full_model_diagnostics\scripts\write_package_manifest.py
```

NEB 的原始有符号端点差与 CatMAP 候选端点最高态有效势垒分别保存；本目录不修改任何科学输入，也不把条件性候选结果升级为正式预测。
