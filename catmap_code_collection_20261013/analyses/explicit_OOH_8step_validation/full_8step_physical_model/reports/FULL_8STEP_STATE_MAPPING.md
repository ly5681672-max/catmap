# 完整显式 OOH 状态映射：候选模型与门控

## 状态映射

原有八步的稳定吸附态保留；完整 NEB 状态不再被压成一个一有效位点物种，而以显式复合态标签表示。当前物种式和组分投影如下，均为根据已有 NEB 状态组成建立的 bookkeeping 映射，**不等同于结构文件已经证明的局部构型或几何占位**：

- NEB1 初态/TS：Ha + OOH + C6H12；NEB1 终态：Ha + C6H12O + O + Hb。
- NEB2 初态/TS：Ha + Hb + O；NEB2 终态：Ha + OH。
- 所以 Ha 在两条路径中均作为显式保留的旁观组分；NEB1 终态生成的 Hb/O 与 NEB2 初态之间有显式状态连接，环氧产物另走吸附态脱附；NEB2 终态的 Ha/OH 再进入原水生成步骤。

为做 CatMAP parser/守恒验证，模型输入对复合态暂用 fragment-count `n_sites`：NEB1 IS/TS=3、FS=4；NEB2 IS/TS=3、FS=2。此数字是可复算的形式占位假设，不是由结构几何验证的物理位点 footprint。因此 `FULL_8STEP_MODEL.mkm` 是 **CONDITIONAL** 候选，不是正式物理模型。

## 计量和库存

CatMAP 0.3.1 对原始八步快照解析 8 个基元步骤、旧九反应候选解析 9 步、当前完整状态候选解析 12 步。候选逐反应原子守恒与形式位点守恒：通过。左零空间维数分别为原始八步 2、旧九反应 1、完整状态候选 2。九反应候选确实少一个库存方向；完整状态方案在代数上恢复了第二个独立不变量。矩阵见 `data/ORIGINAL_8STEP_STOICHIOMETRY.csv`、`data/OLD_NINE_REACTION_STOICHIOMETRY.csv`、`data/FULL_8STEP_STOICHIOMETRY.csv`；相应基见对应 `*_INVENTORY_BASIS.csv`。

映射约定下的库存是：

```text
CatMAP候选矩阵的一组基为 `J1 = θ(Ha*) − θ(Hb*) − θ(OH*) − θ(OOH*)`、`J2 = θ(Hb*) − θ(O*)`。故原研究写法 `I1 = θ(Ha*) − θ(OOH*) − θ(OH*) − θ(O*)` 等于 `J1 + J2`，而 `I2=J2`。对洁净表面零库存类，两者均为0。四个复合态按下列片段投影计算的 `(I1,I2)` 均为 `(0,0)`：NEB1-IS/TS = Ha+OOH+C6H12；NEB1-FS = Ha+Hb+O+C6H12O；NEB2-IS = Ha+Hb+O；NEB2-FS = Ha+OH。该映射保持库存代数，但仍是化学片段记账，不证明真实局域结构或可达性。
```

复合态的 I1/I2 系数按其代表 moieties 投影求和。显式 NEB 反应的这些量保持不变；association/deaggregation 连接也按同一投影保持不变。洁净表面初值对应两常数均为0，但这并不证明模型的稳态根能从该初态动力学可达。

## 能量与动力学边界

Ti–Co 独立态和六个 NEB IS/TS/FS 能量均直接取自隔离的历史 CatMAP 输入快照；没有修改原始 Gibbs。连接能量差由相同能量表的端点求和得到，详见 `data/NEB_CONNECTION_ENERGY_MAP.csv`。现有数值包括 NEB1 `FS→NEB2IS+C6H12O_g = +0.811237 eV`、NEB2 `FS→H2O* = −0.601109 eV`；它们是端点自由能差，不是连接活化能。

NEB1、NEB2 的 TS 保持显式，CatMAP 采用的正/逆能垒分别为 NEB1 `0.556344/0.796030 eV`、NEB2 `0.225994/0.706574 eV`。对未显式 TS 的步骤，本机源码 `catmap/solvers/mean_field_solver.py` 的 `reaction_energy_equations()` 明确令 `G_TS=max(G_IS,G_FS)`；对应正逆端点能垒列在 `data/REACTION_ENERGY_AND_BARRIERS.csv`，故其仍是官方端点近似而非 DFT TS。解析候选总反应自由能为 `−2.088914 eV`，与气相净反应 `H2O2(g)+C6H12(g)→C6H12O(g)+H2O(g)` 的输入能量差闭合。未把端点 ΔG 加到 NEB TS，也未修改 NEB 有符号端点差。

## 结论

现有能量足以建立一个原子守恒、形式位点守恒、并保留两个库存不变量的**条件性拓扑候选**；因此九反应候选不是唯一可表达形式。它还不足以证明正式物理自洽，因为复合态结构投影、Ha 的局域位置、真实 `n_sites`、复合体解离/结合的过渡态及其速率尚未由数据给定。

Ti–Co 已在该候选下尝试零库存稳态，使用 CatMAP 官方速率常数与稳态函数、受库存投影的正对数覆盖度坐标、SciPy 仅作初猜以及 CatMAP `NewtonRoot` 高精度抛光；当前全部有限初猜均未通过残差门槛，无可报告的接受 TOF/覆盖度。尝试记录位于 `runs/tico_provisional_zero_inventory/`。因此，本候选已通过解析和反应化学计量检查，但尚未通过数值稳态门控；更重要的是，复合态几何 footprint 与连接 TS 缺失，不能升级为正式模型。最少缺失数据见 `reports/MINIMUM_MISSING_DFT_DATA.md`。
