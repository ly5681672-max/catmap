# 八步显式 OOH 模型数值门控复核

## 结论先行

旧门控对 Ti、Ti–Hf、Ti–Zr 的“跨精度失败”主要是残差归一化不当造成的误判：跨精度报告把最大物种净残差除以全网络最大净反应速率。近平衡时净速率可能远小于各物种实际正、逆周转通量，这个比值可接近或大于 1，即使每个物种的稳态方程残差相对于其局部正逆总通量已极小。

使用保存根及本机 CatMAP 0.3.1 重新评价而未重新求根后，Ti、Ti–Hf 满足本报告提出的数值通过标准；Ti–Zr 的方程和守恒检查通过，但跨精度 `log10(TOF)` 差为 `2.40095664e-4` decade，超过以 5 ppm TOF 重现精度定义的阈值，因此保留为 `NUMERIC_PARTIAL`。所有数值状态仍是原始八步模型的 `PROVISIONAL` 结果，不代表 NEB 映射已获物理验证。

## 复用数据与执行范围

- 归档根与重放：`H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021\`
- 原质量汇总：`data\MULTISITE_CONDITIONAL_SUMMARY.csv`
- 单位点保存根：`sites\<surface>\constrained_<surface>_result.json`
- 跨精度重放：`sites\<surface>\cross_precision_replay.json`
- 本轮逐物种复核：`C:\Users\Administrator\Documents\codex探索\numerical_gate_audit_20261021\saved_root_local_flux_audit.csv`
- 逐物种复核脚本：`C:\Users\Administrator\Documents\codex探索\evaluate_saved_eight_step_roots.py`

脚本解析归档的 CatMAP 输入，调用 CatMAP 0.3.1 编译后的 `steady_state_function`、速率常数、净反应速率与 TOF API，在已保存覆盖度根上评估；没有调用 `model.run()`、没有重新求根，也没有改动 H 盘项目文件。

## 旧门槛逐项核查

| 旧条件 | 数学/数值含义 | 本次结论 |
|---|---|---|
| 缩放残差 `<1e-50` | 自定义无量纲停止/验收量；对 385 位精度根可达到，但远超一般报告精度所需 | Ti、Hf、Zr 的跨精度根分别为 `1.0947e-56`、`1.0793e-56`、`4.0330e-56`，均已通过；不是三者失败原因。作为统一科学验收阈值过度苛刻且不应强制要求。 |
| 相对残差 `<1e-4` | 旧实现为 `max_i |dθ_i/dt| / max_j |r_net,j|` | 归一化尺度不匹配。它用全网最大净反应速率，而非每个物种的局部正逆总通量；近乎平衡时分母可趋零，导致比值失真。Ti/Hf/Zr 的跨精度值分别为 `5.1704e-4`、`7.2224e-3`、`1.0`，这是三者共同的旧门槛失败项。 |
| 固定绝对残差 `<1e-45` | 有速率单位，绝对量受单位、压力/浓度输入和反应尺度影响 | 三个位点跨精度绝对残差分别为 `2.6635e-56`、`1.8082e-55`、`2.6015e-55`，都通过；不可作为跨模型、跨速率规模的通用门槛。保留作诊断量。 |
| 库存残差 `<1e-300` | 对参数化内生不变量重复作极端精度检验 | I1/I2 由零库存参数化构造满足，且存档时有有限十进制精度；要求 300 位残差没有物理意义。应该核验代数映射及按存档精度计算的舍入界。 |
| 非负覆盖度、位点和为1 | 可行域及位点守恒检查 | 保存根均满足；本轮实算 `|Σθ−1|` 约 `1.2e-81` 至 `4.7e-81`，I1 约 `1e-101` 至 `1e-100`，I2 为0。 |
| 第二精度求解 | 重复性检查 | 旧流程记录了 `Δlog10(TOF)` 和最大覆盖度差，但没有把 TOF/覆盖度一致性作为完整、独立的通过门控。新标准将方程、守恒、覆盖度和 TOF 重现性分别记录。 |

## 三个位点复核证据

| 位点 | 旧跨精度相对残差（旧失败） | 新逐物种最大 `|F_i|/J_i,gross` | `max |F_i|` | 跨精度 Δlog10(TOF) | OOH 覆盖度 | 覆盖度相对误差保守上界 | 新数值状态 |
|---|---:|---:|---:|---:|---:|---:|---|
| Ti | `5.1704e-4` | `5.4313e-56` | `2.6635e-56`（重放） | `1.7178e-21` | `2.0817020e-61` | `<7.7661e-5` | `NUMERIC_PASS` |
| Ti–Hf | `7.2224e-3` | `1.0797e-63` | `1.8082e-55`（重放） | `1.6955e-22` | `1.3796544e-21` | `<2.7495e-46` | `NUMERIC_PASS` |
| Ti–Zr | `1.0` | `1.3769e-48` | `2.6015e-55`（重放） | `2.4009566e-4` | `1.0287055e-21` | `<5.9228e-46` | `NUMERIC_PARTIAL` |

TOF 为保存高精度根重新调用 CatMAP 得到：Ti `5.1513361e-53 s^-1`，Ti–Hf `2.5036393e-53 s^-1`，Ti–Zr `2.5428806e-71 s^-1`。旧相对残差不是这些位点的方程不平衡证据；但 Ti–Zr 的跨精度 TOF 差仍是独立的精度敏感信号，不能抹去。

### 通量接口的重要边界

检查本机 `G:\anaconda\Lib\site-packages\catmap\solvers\mean_field_solver.py`：`get_directional_rates()` 调用 `get_elem_ec()`，后者将速率表达式结果转成 Python `float`。因此不得用该方向接口的两个浮点数相减，来判定 `1e-53` 或 `1e-71` 级的产物净通量；这种相减会丢失有效数字。本次改用 CatMAP 高精度 `get_rxn_rates()` 的净速率，并按官方 `get_turnover_frequency()` 的气相反应符号规则重建 C6H12O 产物速率；三个目标位点重建值与官方 TOF 完全相同（本次精度下相对误差为0）。这证明 API 产物通量的一致性，但不替代跨精度检查，也不证明物理模型正确。

## 全部14位点新旧状态

新规则复核的归档快照为：`data\OLD_VS_NEW_QUALITY_STATUS.csv`。在该快照中，10个位点为 `NUMERIC_PASS`（原8个通过者加 Ti、Ti–Hf），Ti–Zr 为 `NUMERIC_PARTIAL`，Ti–Mo、Ti–W、Ti–Ta 没有已接受根，列为 `NUMERIC_FAIL`。这不是对三个位点的新求解失败声明，而是当前归档中没有通过根可供新门控评估。所有通过只表示数值质量，科学状态仍为 `PROVISIONAL`。

## 输入保护

正式 `.mkm` 与能量输入没有修改。复核前后 SHA-256：

- `alkene_epoxidation.mkm`: `86115CEFB7784FE28AB519271CDFAF0E308B6F6E3A8308B1E90E4318CB497CFE`
- `energies_dft_neb.txt`: `55AD4E94D9C603B6705EC8C92630C9A1CEFF5859CF69865145152375E203CA30`

