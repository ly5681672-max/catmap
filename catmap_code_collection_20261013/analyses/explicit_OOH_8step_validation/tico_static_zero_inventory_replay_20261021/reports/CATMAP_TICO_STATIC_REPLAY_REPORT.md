# Ti–Co 原始八步模型：static 储库零库存重放与条件性 DRC

## 结论

在原始八步 `.mkm` 的 `pressure_mode='static'` 下，CatMAP 0.3.1 重放通过零库存稳态数值门控；与既有浓度探针根的 9 个稳定覆盖度及产物 TOF 完全一致（浮点/存档精度下最大覆盖度差为 0，`log10(TOF)` 差为 0）。基于同一零库存流形、从接受根连续追踪的能量中心差分 DRC 共 44/44 个扰动求解通过数值门控。

**所有新数值仍为 `PROVISIONAL` 条件性结果**：NEB 完整共吸附态与八步网络独立态之间的物理状态映射仍未验证；局部零库存稳态也不等同于已证明从洁净表面初态全局可达。此处没有修改反应网络、能量、位点数或温度。

## 实际输入及储库核对

- 基准反应：原始八步 `alkene_epoxidation.mkm`，解析为 8 条基元反应、9 个稳定吸附物种，CatMAP 0.3.1 `ThermodynamicScaler` + `SteadyStateSolver`。
- 独立运行设置仅把模型运行范围收窄为 `surface_names=['tico']`、温度描述符、`resolution=1`，避免执行原文件的二维描述符网格；原八步反应表达式、`n_sites`、原始能量文件、温度及 static 储库值均未改变。
- CatMAP 实际解析气体顺序下的 static 储库向量为 `[C6H12O_g, C6H12_g, H2O2_g, H2O_g] = [1e-20, 1.0, 0.185, 0.815]`。历史零库存根使用 concentration 探针，但其向 CatMAP 传入的有效向量数值完全相同；本次 static 重放结果与历史根逐项一致。
- 气体数值在 `.mkm` 中以 `pressure` 字段给定，并由 `static_pressure()` 固定读取。这里保留这些输入数值，不将它们重新解释成 bar 或 mol/L。

## 零库存稳态

原八步模型的稳定物种顺序及覆盖度见 `data/constrained_tico_result.json`。关键值：

| 内部物种标识 | 覆盖度 |
|---|---:|
| `C6H12_s` | 0.9895345222678479524 |
| 空位 `*` | 0.01046547773215203653 |
| `OOH_s` | 3.1966965438196039×10⁻²⁰ |
| `Ha_s` | 3.1989328820062164×10⁻²⁰ |
| `O_s` / `Hb_s` | 2.2321519458350×10⁻²³ |

- `I1 = θ(Ha*)−θ(OOH*)−θ(OH*)−θ(O*) = 2.26×10⁻⁵⁷¹`；`I2 = θ(Hb*)−θ(O*) = 0`。
- 位点总和 `=1`；覆盖度非负。
- 8 步净反应速率彼此一致，均约 `1.3236777105111×10⁻³⁶`；CatMAP 的 `C6H12O_g` 净 TOF 为 `1.3236777105111254×10⁻³⁶`，`log10(TOF)=-35.87821774419`。
- 最大独立表面残差约 `1.18×10⁻¹⁰⁰`；相对最大反应通量残差 `8.95×10⁻⁶⁵`。385 位工作精度的独立微扰重放：最大覆盖度差 `2.39×10⁻⁷⁶`，`log10(TOF)` 差 `5.16×10⁻⁴⁹`。

## 两个 NEB 步骤：原始 NEB 与原八步 CatMAP 端点近似分开报告

CatMAP 0.3.1 `templates.rate_constants` 对每个反应使用 `G_TS,eff=max(G_IS,G_FS,G_TS)`，再分别计算 `G_TS,eff−G_IS` 与 `G_TS,eff−G_FS`。所以必须区分全共吸附 NEB 路径与当前输入中独立态所定义的网络端点。

| 路径 | 全共吸附 NEB：ΔG；TS−IS；TS−FS (eV) | 原八步独立态网络 ΔG (eV) | 原八步 CatMAP 实际有效正/逆势垒 (eV) | 解释 |
|---|---:|---:|---:|---|
| 第4步，氧转移 | −0.239686；+0.556344；+0.796030 | +1.8006475 | +1.8006475 / 0 | TS 在该独立态反应的 FS 以下，官方端点最高态规则由 FS 控制；这个有效势垒不是全 NEB 的 +0.556344 eV。 |
| 第6步，Hb/O 氢转移 | −0.480580；+0.225994；+0.706574 | −2.0620670 | +0.225994 / +2.288061 | 该映射下 TS 高于独立态两端，正向值恰与 NEB TS−IS 相同；逆向有效势垒不等于全 NEB 的 +0.706574 eV。 |

因此本次只核对了 CatMAP 实际计算所用的势垒恒等式，没有把端点差异作为 TS 修正，也没有声称映射问题已解决。原八步计算仍保留其输入定义；物理预测效力受映射差异限制。

## 条件性同库存 DRC

官方 `MeanFieldSolver.get_rate_control()` 会通过常规 `get_turnover_frequency()` 重新走未约束覆盖度求解，无法保证保留 `I1=I2=0` 和当前分支。故此处没有把 stock DRC 的未约束结果当成可信值。

下面的 DRC 使用 CatMAP 原生速率常数、方向反应速率和稳态函数，以 `DRC_i=−kBT·∂ln(TOF)/∂G_i` 做 ±能量中心差分；每次扰动均从已接受根做延续，重新在同一精确库存图上求根。步长 0.005/0.010 eV，所有 44 个正/负扰动均通过残差、库存、正覆盖度和位点守恒检查。具体结果见 `data/conditional_drc.csv`；两个步长的差异见同表。结果仍是当前数值根的局部条件性 DRC，不升级为正式机理结论。

幅度较大的局部响应包括 `C6H12_s≈−2.968`、产物态约 `+0.9998`、`Hb_s` 与 `O_s` 约 `+0.9996`、`Ha_s≈+0.5001`、`OOH_s≈−0.4997`。该列表不能单独证明机理速率控制，因为状态映射仍未物理闭合；`OOH-C6H12_s` 扰动 DRC 为 0，与第4步在当前独立态映射下 TS 能量低于最终端点、CatMAP 由终态决定有效势垒相一致。

## 可复现文件与来源哈希

- 执行脚本：`scripts/replay_tico_zero_inventory_static.py`
- 实际运行 `.mkm` 快照：`inputs/tico_probe.mkm`
- 稳态、覆盖度、8 步净通量、残差、库存和详细平衡信息：`data/constrained_tico_result.json`
- 独立精度重放：`data/cross_precision_replay.json`
- 初值/求根尝试：`data/constrained_solver_attempts.json`
- DRC 与 44 条扰动记录：`data/conditional_drc.csv`、`data/conditional_drc_attempts.json`
- DRC 方法边界：`reports/conditional_drc_method.txt`

原始输入 SHA256：

- `alkene_epoxidation.mkm`: `86115CEFB7784FE28AB519271CDFAF0E308B6F6E3A8308B1E90E4318CB497CFE`
- `energies_dft_neb.txt`: `55AD4E94D9C603B6705EC8C92630C9A1CEFF5859CF69865145152375E203CA30`

本次未运行 14 位点批算、未改原始文件、未提交 Git。
