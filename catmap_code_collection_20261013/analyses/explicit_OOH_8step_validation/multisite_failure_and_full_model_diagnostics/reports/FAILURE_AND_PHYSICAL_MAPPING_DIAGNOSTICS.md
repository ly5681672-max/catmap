# 八步失败位点与完整模型物理映射诊断

## 摘要

- 批次内候选模型：CatMAP 官方解析为 9 个基元反应；原始显式网络解析为 8 个基元反应、八个化学阶段。
- 表面计量矩阵的守恒零空间维数：原始网络 2，九反应候选 1。原始网络的两个库存不变量与候选网络的单个 `Ha−OOH` 不变量不等价；因此候选网络不是已证明的严格动力学粗粒化。
- 原始 TS 低于至少一个端点的 NEB 路径：7 条，涉及 7 个位点。每条记录均保留原始 TS−IS、TS−FS 和 CatMAP 端点最高态正逆势垒。
- 10 个已收敛候选中主导状态分布：adsorbed 1-hexene-labelled state: 3; NEB1 final composite state: 6; adsorbed water-labelled state: 1。

## 失败位点诊断与恢复尝试

`failed_site_root_cause_diagnostics.csv` 保留能量输入跨度、已有尝试残差、求解异常证据及可证范围。结论仅限数值：

- Ti–Zr：既有批次在求解诊断前中断；本轮恢复流程只完成 `uniform_0` 一个初猜，尚未运行 donor 初猜，也未完成完整验收；不能把它归类为已证明的 CatMAP 数值失败根因。
- Ti–Ti、Ti–Ni：多个初猜在 120/240 位精度下有很大的相对表面残差，属于未找到稳态根；当前证据不证明不存在根。变量仍采用正覆盖度 log 参数化、形式位点守恒和由计量矩阵推导的库存约束。
- Ti–W：尝试日志出现数值奇异矩阵且 240 位任务超时。其 NEB1 原始 TS−IS = −0.566651 eV；端点最高态法因此给候选正向有效势垒 0 eV。该物理异常标记与求解奇异是两类问题，不能互相替代解释。

本轮 Ti–Zr donor-seed 恢复脚本首个 uniform 初猜达到缩放残差 `4.5881e-86`，但旧相对残差门槛仍为 1，未通过。该恢复进程达到本轮时间预算后停止；没有生成可接受根或 TOF。SCIPY 初始 Jacobian 条件数约 `3.58e5`，该数值只针对这个初猜，不能外推为全局条件数。日志/脚本保留了该次有限尝试；未降低门槛。

## NEB 异常与势垒语义

7 个要求核对的位点均由自动规则得到明确标记：Ti 和 Ti–Zr 为 `TS_BELOW_FS`；Ti–W、Ti–Ta、Ti–Mo、Ti–V、Ti–Cr 为 `TS_BELOW_IS`。这是原始有符号比较，不将负值钳制或重写为真实零势垒。CatMAP 候选有效势垒按 `max(G_IS,G_TS,G_FS)` 计算，分别列出正、逆方向；端点值不是新 DFT TS。

## 主导覆盖度与候选机制解释

候选中 NEB1 终态复合态覆盖接近 1 的位点确实包括 Re、Mo、Mn、Cr、V、Ta；Hf 的最大覆盖是水标记态；Ti、Co、Fe 的最大覆盖是烯烃标记态。此分布与各自能量表和候选网络平衡相容，但复合态是粗粒化的一有效位点表达，不证明真实几何占位。

| 位点 | 候选模型最大覆盖状态 | 最大覆盖度 | OOH* 覆盖度 | 空位覆盖度 | 产物 TOF (s⁻¹) |
|---|---|---:|---:|---:|---:|
| ti | adsorbed 1-hexene-labelled state | 0.968735 | 1.94148e-33 | 0.0311877 | 2.53476e-53 |
| timn | NEB1 final composite state | 0.996449 | 8.15944e-12 | 5.1401e-05 | 1.50241e-12 |
| tihf | adsorbed water-labelled state | 0.643188 | 1.20505e-21 | 0.0141607 | 3.00439e-30 |
| tire | NEB1 final composite state | 0.982759 | 5.16222e-09 | 7.2021e-05 | 8.21174e-07 |
| timo | NEB1 final composite state | 0.985357 | 5.10132e-10 | 9.45694e-07 | 2.75099e-10 |
| tiv | NEB1 final composite state | 0.989618 | 1.03158e-13 | 0.000349246 | 7.36301e-16 |
| tita | NEB1 final composite state | 1 | 1.28982e-09 | 7.94808e-13 | 7.28314e-24 |
| ticr | NEB1 final composite state | 0.999689 | 4.1921e-12 | 2.66074e-07 | 1.15028e-14 |
| tico | adsorbed 1-hexene-labelled state | 0.989535 | 2.26714e-20 | 0.0104655 | 3.53064e-27 |
| tife | adsorbed 1-hexene-labelled state | 0.97985 | 3.27929e-22 | 0.00570208 | 7.31453e-31 |

Re 与 Mo 的 NEB1 终态高覆盖同时伴随较高候选 TOF；所以“NEB1终态覆盖高必然堵塞空位并降低 TOF”不成立于这组模型数据。Ta 的空位覆盖极低，可能在候选模型内限制需要空位的连接，但仅凭覆盖度不能定为速率控制。跨位点 TOF 差异同时依赖各能量输入、端点最高态势垒、气相储库和候选连接步骤。

**连接动力学敏感性边界：**现存 `MULTISITE_LINK_ENDPOINT_SENSITIVITY.csv` 中，Ti–Co 的连接端点/假设敏感性来自历史条件性计算，Ti–Fe 仅有少量复用的端点能量敏感性；Re、Mo 对应记录均为 `NOT_RUN`。本轮没有为 Re/Mo 引入任意连接 TS 或借用 Ti–Co 势垒，因此不能据现有数据断言 Re/Mo 的 TOF 排序对连接动力学稳健。表中 `conditional_sensitivity` 也不是 CatMAP 官方 DRC。

## 原始八步与九反应候选的数学关系

实际 CatMAP parser 解析所得计量矩阵和 nullspace 系数见 `data/`。原始网络库存基底等价于：

```text
I1 = θ(Ha*) − θ(OOH*) − θ(OH*) − θ(O*)
I2 = θ(Hb*) − θ(O*)
```

九反应候选的库存只留下：

```text
I' = θ(Ha*) − θ(OOH*)
```

候选将 Hb*、O*、OH*、独立环氧吸附态及其库存变化并入 NEB 复合态/连接，不能从原始网络严格推导出此减少。两网络按所声明分子式和每个吸附物种 `n_sites=1` 做逐反应原子与形式位点计数均通过；这只验证代数平衡，不验证复合态结构身份、真实几何占位、能量连续性或连接路径的 DFT 势垒。完整矩阵、逐反应检查与状态映射见 `data/original_8step_surface_stoichiometric_matrix.csv`、`data/candidate_9reaction_surface_stoichiometric_matrix.csv`、`data/network_balance_checks.csv`、`data/original_to_candidate_reaction_mapping.csv`。

## 结论与门控

现有数据可以构成原子/形式位点守恒的九反应**条件性候选**，但没有足够依据将它认定为原始八步网络的严格粗粒化。NEB1/NEB2 原始三态可在候选中保留；缺失的连接 TS、Hb/O/OH 独立库存、复合态真实 `n_sites` 以及从洁净表面初态的可达性仍阻止物理升级。现在不能据此开展正式 14 位点 TOF/DRC 或正式构效拟合。

`MULTISITE_CALCULATION_STATUS.csv` 的 `reused_existing_result` 对部分新算位点为 True，确有口径歧义。旧表未覆盖；本次数据整理应区分“本批 CatMAP 重求解”“复用历史稳态”和“仅复用输入/预处理产物”，对应于 `calculation_provenance_audit.csv`。

本次另生成 `data/multisite_candidate_master.csv`，按 `surface_name` 将14个目标位点状态、已有覆盖度/TOF、NEB1/NEB2有符号端点差与候选有效势垒、旧的条件性TS敏感性及连接端点敏感性状态合并。空字段表示来源结果未提供该项，不以0替代；该表仍属于候选模型数据。
