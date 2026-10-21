# Ti–Co 显式 OOH 候选模型：速率机制与条件性敏感性

**状态：PROVISIONAL。** 本目录只诊断既有 Ti–Co 候选，不是物理映射已验收的正式八步预测。所有输入/情景都隔离保存；未修改原始八步 `.mkm`、Excel、DFT/NEB、CatMAP 安装文件或历史结果。

## 本轮结论

- 333.15 K 条件下严格零库存约束候选根复现：产物速率约 `3.53063694010×10^-27`（CatMAP 本模型速率单位）；覆盖度位点和为 1，且模型明确施加 `θ(Ha*)=θ(OOH*)`。
- 稳态正逆通量显示，过氧化氢活化与“吸附氢/过氧基/烯烃到完整 NEB1 初态”的连接都与同一串联通量相连。当前均场公式使第4反应正向通量正比于 `θ(Ha*)θ(OOH*)θ(C6H12*)`，在模型约束下为平方形式；这只是模型假设，不是实验反应级数。
- 按 ±0.005/±0.010 eV 重新求稳态的 NEB-TS 条件性 DRC 极小；反而两个连接步骤的假设性附加势垒敏感性约 0.52。后者并非 DFT DRC，且强烈依赖缺失 TS 的情景参数。
- CatMAP 官方 ODE 从空表面在 300 步后未达稳态；这不证明物理不可达。数值根本身也不证明由洁净表面动力学轨迹可达。

## 文件导航

- `FORWARD_REVERSE_RATE_AUDIT.csv`：9步 `kf/kr`、正逆通量、高精度净通量、端点自由能/有效势垒与稳态覆盖度；保留 CatMAP 官方方向率 API 的双精度输出及其精度限制说明。
- `TICO_PROVISIONAL_DRC.csv`：NEB1/NEB2 显式 TS 的中心差分条件性 DRC。
- `EXPLICIT_TS_PERTURBATION_CASES.csv`：8 个扰动点的输入/稳态/日志索引。
- `ENDPOINT_BARRIER_CONTROL.csv`：7个无显式 TS 连接的 δ=0、0.005、0.010 eV 情景及单侧敏感性。
- `HA_OOH_MEAN_FIELD_ASSESSMENT.md`、`DRC_NUMERICAL_VALIDATION.md`、`CLEAN_SURFACE_REACHABILITY_DIAGNOSIS.md`、`MULTI_CATALYST_8STEP_ROADMAP.md`：机理边界、敏感性质量、空表面可达性和扩展路线。
- `runs/`：各条件 CatMAP 输入快照、收敛根、覆盖度、速率、尝试记录和日志。
- `candidate_tico.mkm`、`candidate_energies.txt`、`candidate_reaction_energy_barriers.csv`：本分析自包含候选输入快照；不得替代正式模型。
- `INPUT_SHA256.txt`、`OUTPUT_SHA256.txt`：来源输入与本目录结果清单（结果清单覆盖本包除两份清单外的文件）。
- `validation_summary.json`：自动质量门控摘要。

## 复现

在本机 Python/CatMAP 0.3.1 环境运行：

```powershell
$env:PYTHONIOENCODING='utf-8'
G:\anaconda\python.exe .\run_conditional_sensitivity.py --base-run <候选基线运行目录> --solver-script .\run_candidate_zero_inventory_strict.py --output .
G:\anaconda\python.exe .\run_rate_audit.py --model .\runs\strict_baseline_replay\candidate_tico.mkm --coverages .\runs\strict_baseline_replay\candidate_steady_state_coverages.csv --output . --barriers .\candidate_reaction_energy_barriers.csv --precision 220
G:\anaconda\python.exe .\finalize_sensitivity_tables.py
G:\anaconda\python.exe .\validate_kinetic_control.py
```

当前验证摘要 `validation_summary.json`：23/23 条件性稳态求解通过本目录严格数值门控；方向通量、库存和位点约束核验通过。这个数值 PASS **不提升物理状态为正式**。
