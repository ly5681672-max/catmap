# 显式 OOH 多催化剂候选工作流：首个非 Ti–Co 数值验证

## 结论摘要

- 成功复用归档的 Ti–Co 候选基线；未重跑其既有 23 组扰动。
- 依据项目自己的 Excel 转换脚本和 DFT/NEB 总 Gibbs 表，构建并实际运行 Ti–Fe 条件性候选。CatMAP 0.3.1、220 位工作精度下数值门控通过，产品 `C6H12O_g` TOF 为 `7.3145289447224763e-31 s^-1`。
- 该结果不是正式八步物理预测：完整 NEB 复合态的单有效位点粗粒化、独立吸附物种到共吸附 NEB 初态的映射，以及缺失连接 TS 的端点势垒仍未被 DFT 结构/动力学证据验证。
- 本机 CatMAP 0.3.1 将输入式中的 `->` 与 `<->` 统一解析成初态/过渡态/终态三组，不保留箭头方向；其速率常数函数对全部 elementary reaction 计算正、逆常数。因此候选第2步的单箭头并未实现单向动力学。现有 Tico 历史包中“单向第2步”解释需后续单独校正，历史文件本轮未改。

## 已复用及源数据

1. Ti–Co 复用文件：`analyses/explicit_OOH_8step_validation/candidate_tico_kinetic_control/`。其 `validation_summary.json` 记载 23/23 数值案例通过，最大相对表面残差 `4.81e-62`，最大位点/Ha−OOH误差为0；其结果仍标为旧候选 `PROVISIONAL`。本轮不重跑这些案例。
2. Ti–Fe 输入由 `H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\DFT_formation_catmap.txt`、`DFT_TS_relative.csv`、该项目的 `己烯环氧化.xlsx` 和 `convert_dft_to_catmap.py` 生成；来源快照及 SHA256 位于 `inputs/source_snapshot/`。权威项目工作簿 SHA256：`90f1700457ece5a53335ef771e88e24134c00afab57b3d8ec15a18bb307f38bd`。
3. Ti–Fe 使用工作簿 `catmap数据集` 第51行同位点 clean-slab G `−870.962186 eV` 及既有元素参考，将六个 NEB 总 Gibbs 换成候选 formation energy；该代数转换可复算，但不证明结构、超胞、真实几何占位与独立吸附态完全等价。

### Ti–Fe 候选 NEB 输入与原始有符号端点差

| 状态 | 候选 formation energy (eV) |
|---|---:|
| NEB1 IS / TS / FS | −33.597103 / −33.437687 / −35.776086 |
| NEB2 IS / TS / FS | −7.029853 / −6.8568715 / −7.646919 |

对应原始 signed `TS−IS`/`TS−FS` 分别为 NEB1 `+0.159416 / +2.338399 eV`、NEB2 `+0.1729815 / +0.7900475 eV`。未截断、未人为改势垒。NEB1复合态组成按 `C6H14O2`，NEB2按 `H2O` 计入元素参考；真实结构相同/占位对应关系仍标记 UNVERIFIED。

## 网络映射与守恒

条件性候选用九条 CatMAP elementary reaction 表示原八个化学阶段：四个独立前驱态先汇入 NEB1 完整 IS，之后为 NEB1 IS–TS–FS；NEB1 FS 接到 NEB2 IS 并释放环氧产物；随后为 NEB2 IS–TS–FS、接入独立吸附水态并脱附。完整反应式、状态映射和边界见 `CANDIDATE_NETWORK_MAPPING.md`。

CatMAP 实际解析出的九个稳定吸附物种为 H2O2、C6H12、Ha、OOH、H2O、NEB1IS/FS、NEB2IS/FS；NEB1TS/NEB2TS 不作为稳态覆盖度。候选中所有这些物种均设 `n_sites=1`，因此形式位点平衡为

`theta_* + theta_H2O2 + theta_C6H12 + theta_Ha + theta_OOH + theta_H2O + theta_NEB1IS + theta_NEB1FS + theta_NEB2IS + theta_NEB2FS = 1`.

由候选计量矩阵精确求得其零库存不变量为 `theta_Ha − theta_OOH = 0`。它不是原始八步完整网络的库存不变量。完整八步网络有两个不变量，且显式含 Hb/O/OH 状态；当前九步模型没有证明严格动力学粗粒化，也不能把该等式外推到原始八步模型。

## Ti–Fe 数值结果（候选条件）

| 指标 | 实测数值 / 状态 |
|---|---:|
| `theta_*` | 0.0057020816259998064 |
| `theta_C6H12*` | 0.9798503514592182 |
| `theta_OOH* = theta_Ha*` | 3.2792896617874005×10^-22 |
| `theta_NEB2FS*` | 0.01439676820784059 |
| `theta_H2O*` | 5.079870027960455×10^-5 |
| `TOF(C6H12O_g)` | 7.3145289447224763×10^-31 s^-1 |
| 接受的相同分支根 | 2；覆盖度最大绝对离散约 8.21×10^-95 |
| 最大相对表面稳态残差 | 3.02×10^-54 |
| 位点平衡误差 / 库存不变量误差 | 4.43×10^-221 / 0 |
| 九步净通量与产物 TOF 最大相对差 | 3.02×10^-54 |

所有覆盖度非负，形式位点总和为1。上述是候选模型数值内部一致性，不证明物理模型闭合或从洁净表面可达。

## 详细平衡与第2步方向语义

使用本机 CatMAP `solver._kB=8.617332478e-5 eV K^-1`，逐步由解析形成能计算 `DeltaG=G_products−G_reactants`，并验证 `kBT ln(kf/kr)+DeltaG`。九步最大偏差低于 `1e-8 eV`，`tife_detailed_balance_audit.csv` 全部为 CHECKED。第2步输入虽然写 `H2O2* + * -> Ha* + OOH*`，实际 reverse rate constant 非零；因此这组稳态并非研究者拟定的“单向第2步”动力学。没有为实现单向而私自置零逆速率。

## 条件性 DRC 与连接敏感性

每个扰动重新用 CatMAP 能量/速率函数和相同 `Ha−OOH=0` 约束求稳态，δ=±0.005 eV，分支对照均通过。因为本机官方 `get_rate_control` 对能量扰动后会自行重新求覆盖度，不能证明维持本候选的零库存子空间，本轮不把其未经约束的内置输出冒充该分支 DRC；使用约束重求的有限差分并明确标为条件性 DRC：

- NEB1 TS：`6.85×10^-41`；NEB2 TS：`1.10×10^-40`，均数值上近零。
- NEB1 FS、NEB2 FS 为无显式 TS 连接端点的能量敏感性，不是 DRC；NEB2 FS 的局部斜率约 `1.262 eV^-1`。它不等于缺失连接势垒敏感性，也不能给出真实连接 TS。

## 多位点数据完整性和范围

`data/site_readiness.csv` 对原声明目标 14 个表面（Ti + 13 个 Ti–M）逐个检查。14/14 当前表格具有六个独立输入字段和 NEB1/NEB2 六个 IS/TS/FS 数值；这仅为“数据行存在”，不代表能量映射、TS 物理有效性、共同结构参考已通过，也不代表运行完成。Ti–Nb 在源表中额外存在，仍作为额外候选，不混入原目标清单。Ti–Zr未运行、未解决。

源表中的 signed NEB TS−IS 有若干负值：Ti–W NEB1 `−0.566651 eV`、Ti–Ta NEB1 `−1.174636 eV`、Ti–Mo NEB2 `−0.012820 eV`、Ti–V NEB2 `−0.024677 eV`、Ti–Cr NEB2 `−0.1931555 eV`。这些记录原样保留；不能把它们自动解释为有效 TS，也不能静默截成零。相较之下，所选 Ti–Fe 两个 TS 均高于端点，适合作为数值工作流首测，但仍有结构映射限制。

## 正式计算就绪判断

**尚未达到正式八步 TOF/DRC 与 OOH–TOF 构效拟合条件。** 最小阻碍不是缺数值行，而是：

1. 独立吸附态到含 Ha 旁观物的完整 NEB1/NEB2 IS/FS 的结构、超胞及位点计量映射没有逐位点证明；单有效位点 NEB 状态是明确粗粒化假设。
2. 候选前驱物→NEB1IS、NEB1FS→NEB2IS+产物、NEB2FS→H2O* 等无显式 TS 连接的真实动力学参数缺失；端点势垒近似只支持条件性探索，连接自由能不是活化能。
3. CatMAP 箭头语义需单独解决：当前实际模型中第2步不是单向；若科学模型必须单向，需通过官方 0.3.1 机制/API 或明确动力学实现证明如何设置，不能靠改速率常数数组。
4. 原八步显式网络的库存类 `I1/I2` 与当前粗粒化候选 `Ha−OOH=0` 不等价；没有完整状态的可达性证明。

故当前可以继续按位点复制候选**数值筛查**并逐个标 `PROVISIONAL`，但不应将其混入正式多催化剂 OOH–TOF/DRC/SISSO 数据集。下一里程碑是确定能够在不丢失原八步库存的情况下纳入完整共吸附 NEB 态的 CatMAP 网络/占位表达，并为缺失连接找到真实动力学证据或由研究者批准的清晰条件参数范围。未开展全位点批算。
