# explicit_OOH_8step 多催化剂候选工作流

此目录记录从既有 Ti–Co 候选复用流程并以 Ti–Fe 做首个非 Ti–Co 验证的结果。所有新计算仅为 `PROVISIONAL`；本目录没有改动原始 Excel、DFT/NEB、正式模型、CatMAP 安装源码或历史候选结果。

## 目录

- `inputs/source_snapshot/`：用于 Ti–Fe 转换的源文件快照与 SHA256。
- `runs/tife_candidate/`：Ti–Fe 单个位点候选模型、输入、220 位精度根结果、覆盖度、九步速率和能量映射。
- `runs/tife_sensitivity/`：NEB1/NEB2 TS 及 NEB1FS/NEB2FS 端点 ±0.005 eV 的独立 CatMAP 重求案例及日志。
- `data/`：14 位点数据完整性表、Ti–Co/Ti–Fe 可追溯候选数据表、详细平衡审计、条件性敏感性表。
- `reports/`：原八步到九反应候选映射与当前数值/物理门控报告。
- `scripts/`：从只读源表准备模型、运行 CatMAP、做扰动重求、生成数据清单和审计。
- `tests/`：对输出位点平衡、库存、通量、详细平衡和数据清单的轻量回归测试。

## 复算顺序

在项目根目录以本机 Python/CatMAP 0.3.1 运行：

```powershell
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\prepare_site_candidate.py --surface tife
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\run_site_candidate.py --model .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\runs\tife_candidate\candidate_tife.mkm --precision 220
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\audit_tife_rate_thermodynamics.py
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\run_tife_conditional_sensitivity.py --precision 220
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\build_site_readiness.py
G:\anaconda\python.exe .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\scripts\assemble_provisional_dataset.py
G:\anaconda\python.exe -m unittest discover -s .\analyses\explicit_OOH_8step_validation\multisite_candidate_workflow\tests -v
```

准备脚本读取项目转换脚本指定的 `H:\catmap\catmap\lyq_alkene_epoxidation\己烯环氧化.xlsx`。不同路径的同名 Excel 不可混用。若输入源哈希与既有快照不同，脚本拒绝覆盖不一致快照。

## 重要口径

- 正式主线仍是原八步显式 OOH 网络。此处九反应候选将两段完整共吸附 NEB 放进表面状态序列，并把若干缺失 TS 的连接交给 CatMAP 端点近似；它不是原八步的已证严格粗粒化。
- 稳态零库存条件从实际候选化学计量矩阵左零空间推导；当前条件得到 `theta_Ha−theta_OOH=0`。不要手动套用到原始八步的两个不变量或其他网络。
- CatMAP 输出 `NUMERICAL_PASS_PROVISIONAL_PHYSICS` 表示数值根及通量检查通过，不表示物理验证通过。
- CatMAP 0.3.1 本机解析会把单箭头/双箭头都压成状态序列；检查发现第2步仍有逆速率常数。不要把 `.mkm` 中 `->` 文字当作已落实的不可逆动力学。
- 显式 NEB TS 的条件性 DRC 是在同一约束子空间分别重求后中心差分；无显式 TS 的端点能量敏感性不等于真实缺失势垒 DRC。
- `site_readiness.csv` 只判断输入行是否存在；不表示批量计算或物理门控已通过。TiNb 在数据表存在但不在原14目标列表；TiZr本轮未运行。

## 本轮状态

已实际运行 CatMAP 0.3.1 Ti–Fe 候选稳态与 8 个敏感性变体，详细平衡审计 9/9 通过；5 项输出测试通过。正式八步物理映射、原始初态库存类可达性、完整目标位点计算及正式 DRC 尚未通过或未执行。详见 `reports/MULTISITE_CANDIDATE_RESULTS.md`。
