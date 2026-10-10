# CatMAP 333.15 K 单温度多初值恢复

## 研究条件保持不变

恢复程序只在 **333.15 K、pressure=1.0** 求解。以下内容保持原样：

- `energies_dft_neb.txt`
- `DFT_TS_relative.csv`
- `alkene_epoxidation.mkm`
- 有效气相活度和八步反应机理

550 K、750 K 曾经用于数值桥接，但它们不是研究条件，也不再由恢复程序生成或读取。旧的 `analysis_ooh/recovery_v2/` 目录保留为历史证据，不参与新的正式恢复。

## 单温度多初值方案

对每个未通过验证的表面，在同一个 `[333.15 K, 1.0]` 描述符点生成以下初始状态：

1. 空位主导；
2. OOH* 主导；
3. OH* 主导；
4. O* 主导；
5. C6H12* 主导；
6. H2O*、Ha*、Hb* 主导的合理替代初值。

CatMAP 默认 `numbers_solver` 使用平方 numbers 参数。程序将覆盖度初值转换为对应的 `numbers_map`，每个初值独立运行，不传播前一个初值的目标解。

每个初值进行 260/360/460 位精度的连续复算。只有同时满足以下条件才算有效：

- 八步基元通量闭合；
- 气相物料守恒；
- 覆盖度归一化且非负；
- 稳态残差相对于最大通量达标；
- 连续两级精度的 `log10(TOF)` 差值不超过 0.02；
- 至少两个不同初始状态得到一致的有效根。

不同初值产生不同有效根时，标记为 `different_valid_steady_state_branches`，不自动选择其中一支。

## Windows 运行

```powershell
cd H:\catmap\catmap
git -c http.proxy=http://127.0.0.1:6696 -c https.proxy=http://127.0.0.1:6696 pull --ff-only origin main
cd lyq_alkene_epoxidation\2-creating_microkinetic_model
$env:PYTHONUTF8='1'
python -m unittest -v test_ooh_numerics test_ooh_recovery
python run_ooh_analysis.py --recover --recover-surfaces titi timn ticr
python run_ooh_analysis.py --recover
```

`PYTHONUTF8=1` 只用于避免 Windows 默认代码页遮蔽 CatMAP 的中文错误日志；它不改变反应、DFT 能量或求解器判据。

## 输出

- `analysis_ooh/recovery_single_temperature/seed_definitions.csv`：初始状态、物种顺序和覆盖度向量；
- `analysis_ooh/recovery_single_temperature/seed_trials.csv`：每个初值的精度复算；
- `analysis_ooh/recovery_single_temperature/recovery_results.csv`：跨初值最终判定；
- `analysis_ooh/recovery_single_temperature/target/`：每个初值和精度级别的 CatMAP 日志、pkl 和结果。

只有通过全部质量门槛的表面才会更新 `dft_baseline.csv`、`OOH_TOF_results.csv` 和拟合图。代码运行成功不等于科学结果验证成功。
