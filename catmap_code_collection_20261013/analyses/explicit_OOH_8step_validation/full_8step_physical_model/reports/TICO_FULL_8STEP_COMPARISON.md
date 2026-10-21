# Ti–Co：原始八步、九反应候选与完整状态候选对照

## 结论

现有输入能构建一个保留全部原始稳定吸附态、显式保留 OOH*、保留 NEB1/NEB2 真实过渡态，并在**片段记账假设**下恢复两个零库存不变量的 CatMAP 候选网络。它可被 CatMAP 0.3.1 解析，12 个基元步骤均通过原子和形式位点守恒。它尚不是正式物理模型：共吸附结构的真实几何位点数、状态构型等价性和多个连接过程的真实活化能仍未验证；一次零库存条件性稳态求解未通过数值门槛。

## 三种网络的区别

| 网络 | 状态表示 | 库存不变量 | 物理边界 |
|---|---|---|---|
| 原始八步 | H2O2*、Ha*、Hb*、OOH*、C6H12*、C6H12O*、O*、OH*、H2O* 等稳定态；NEB1、NEB2各一个显式TS | 两个独立库存 | NEB完整共吸附态端点能量与原网络的独立吸附态组合能量不一致，不能直接把其TS−IS移植到独立端点反应 |
| 旧九反应候选 | 通过合并/删除部分 O/H 稳定态简化网络 | 仅一个已施加约束的库存关系 | 未证明严格动力学粗粒化；不能作为完整八步等价模型 |
| 本次完整状态候选 | 保留原有稳定态；新增NEB1-IS/TS/FS与NEB2-IS/TS/FS状态标签；连接用端点过程 | 两个不变量都在片段投影下保留 | parser/原子/形式位点通过；复合态位点footprint与连接动力学是条件假设；Ti–Co零库存求根尚未通过 |

## 端点能量与势垒

Ti–Co NEB1：IS/TS/FS = `−33.720385/−33.164041/−33.960071 eV`；`ΔG(FS−IS)=−0.239686 eV`；真实显式TS给出的正/逆势垒为 `0.556344/0.796030 eV`。

Ti–Co NEB2：IS/TS/FS = `−5.481063/−5.255069/−5.961643 eV`；`ΔG(FS−IS)=−0.480580 eV`；正/逆势垒为 `0.225994/0.706574 eV`。

旧独立态反应能与 NEB 端点差分别相差 `+2.0403335 eV`（NEB1）及 `−1.5814870 eV`（NEB2）。这些差值不是可直接加到过渡态上的修正量。本次候选通过显式完整 NEB 状态避免把 NEB 端点强行等同于独立吸附态总和；其余连接自由能由同一张现有输入表的化学计量端点和计算，见 `../data/NEB_CONNECTION_ENERGY_MAP.csv`。

## 求解尝试

运行脚本调用 CatMAP 原生反应解析、Scaler/速率常数和稳态函数；对零库存状态使用正对数覆盖度坐标，SciPy 只产生局部初猜，随后使用 CatMAP `NewtonRoot`。所有12步的 CatMAP 原生 `kf/kr` 与 `exp(−ΔG°/kBT)` 检查通过，最大相对误差 `1.03×10⁻¹³`。旧九反应稳态只作为数值初猜来源，已经按本候选两条库存关系投影，不复用其TOF。4组有限初猜均未满足严格残差门槛，故没有可接受的新候选稳态、TOF或覆盖度。失败详情见 `../runs/tico_provisional_zero_inventory/solver_attempts.csv` 与 `run_summary.json`。

## 可复现产物

- 候选输入：`../inputs/FULL_8STEP_MODEL.mkm`（文件头声明 CONDITIONAL）。
- 逐反应计量与库存：`../data/FULL_8STEP_STOICHIOMETRY.csv`、`../data/FULL_8STEP_INVENTORY_BASIS.csv`。
- 原始八步对应矩阵：`../data/ORIGINAL_8STEP_STOICHIOMETRY.csv`、`../data/ORIGINAL_8STEP_INVENTORY_BASIS.csv`。
- CatMAP 端点/TS能量与势垒：`../data/REACTION_ENERGY_AND_BARRIERS.csv`。
- 原始输入哈希：`../data/SOURCE_HASHES.csv`。
