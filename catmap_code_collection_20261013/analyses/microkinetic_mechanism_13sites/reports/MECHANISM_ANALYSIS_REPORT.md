# 13个Ti–M位点 CatMAP 微观动力学机制结果

## 范围与来源

- 计算：本报告复用已完成的 CatMAP 0.3.1 严格容差六步稳态输出；未重跑稳态基线。
- 位点：13个已验证位点；Ti–Zr未纳入。
- 条件：333.15 K；固定基线储库输入 H2O2=0.185、C6H12=1.0、H2O=0.815、C6H12O=1e-20；模型中的压力描述符=1.0；其单位/物理量遵循基线配置，不在此处重新解释。
- 机理：原六步反应网络不变；第2/4步保留显式NEB TS1/TS2，第1/3/5/6步为 CatMAP 端点势垒近似。
- DRC：5个位点复用既有官方 `rate_control` 输出；另8个位点通过 CatMAP 0.3.1 官方 `ReactionModel.run(recalculate=True)`、仅请求 `rate_control` 的隔离补充运行取得。该补充没有改动能量、机理、温度或储库。
- 主要输入与输出路径及SHA256见 `../data/source_SHA256.csv`；补充DRC运行日志、配置、能量表位于 `../data/official_drc_runs_supplemental/`。

## 数值验收摘要

- 严格稳态：13/13 为 `CONVERGED`。
- 最大覆盖度位点守恒误差：1.000e-28；覆盖度检查同时保留高精度十进制原值。
- 第3步净速率与 C6H12O 气相净TOF最大差：0.000e+00（以导出十进制值比较）。还用CatMAP CSV中每步反应两端的实际物种计量系数求和核算了六步网络产物速率；结果见 `TOF_flux_consistency.csv`。
- 原方向速率CSV表头有3列而数据行实际为2列；脚本按原始数据行结构解析。方向API输出含负值：forward 5 行、reverse 5 行；最大缩放闭合差为 1.179e+00。本机CatMAP 0.3.1 `mean_field_solver.py` 中 `get_empty_site_cvgs()` 将覆盖度转成Python `float` 后计算空位，`get_elem_ec()` 又将速率常数、覆盖度和结果转成`float`；当空位接近零时会产生负空位/负方向输出。因此第1、6步方向API数值在这些体系中不能作为可信单向物理通量；原始值保留并在诊断CSV中标记。此问题不推翻高精度稳态净速率、覆盖度或官方TOF。
- DRC输出覆盖：13/13位点 × 7个表面态/过渡态参数；DRC总和范围=-0.895431 至 0.998818。总和偏离1时不强制归一化，保留CatMAP原值。
- 已有TS1–TS2 50×50 CatMAP缓存的pickle中确实含2500个`coverage_map`节点，且与既有TOF节点按描述符坐标一对一匹配；本报告复用该缓存，没有重跑网格。覆盖度向量按CatMAP日志的5个`adsorbate_names`顺序加末尾空位`s`核验，逐节点总和为1。峰值及TS1/TS2轴向两侧代表节点详见`ts1_ts2_grid_peak_and_sides.csv`；完整覆盖度长表和Origin矩阵一并导出。虚拟网格数据与13个真实催化剂数据分开保存。
- 速率常数：基线CSV没有直接导出每步正/逆速率常数，因此 `reaction_rate_constants_13sites.csv` 将其记为 `NOT_EXPORTED_IN_EXISTING_BASELINE_CACHE`；未由净速率反推或填造。

## 机制结果摘要

- 覆盖度最大稳定吸附态按位点列于 `coverage_dominant_species.csv`（见覆盖平衡表）；代表性最高覆盖的前三个位点：EpoxideHHO*—Ti–W, Ti–Ta, Ti–Mn；HHO*—Ti, Ti–Co, Ti–Fe；OHHa*—Ti–Fe, Ti–Ti, Ti–Ni。请结合逐位点覆盖度，不把最大覆盖物种直接等同于速率控制步骤。
- 已有虚拟网格离散最大点在 TS1=-33.314735 eV、TS2=-7.673499 eV，log10TOF=4.8219；峰处θ(EpoxideHHO,H2O,HHO,OHHa,Rpair,*)=(0.123,0.477,5.82e-05,0.00475,0.395,5.24e-08)。沿TS1和TS2两轴各±2格时，H2O/Rpair/EpoxideHHO覆盖重新分配（见5点及2500节点明细）；最近真实催化剂是Ti-Ni，描述符平面欧氏距离=0.415 eV，真实log10TOF=-1.630。真实位点并未达到虚拟网格峰附近的预测活性，不能把该离散峰等同于实际催化剂的可实现性能。
- 全部13个位点的产物TOF、六步高精度净通量、方向API原始输出、有效正逆势垒和官方DRC分别见数据表与图；方向速率图把负API值遮罩为灰色，不解释为真实通量。
- 有效势垒图中第1/3/5/6步明确标为 `CatMAP_endpoint`，不能解释为DFT独立过渡态能垒；第2/4步为输入的显式NEB TS势垒。
- TS1−Rpair与TS2−HHO的原始有符号差值逐位点见 `explicit_TS_signed_energy_differences.csv`；负值原样保留，和CatMAP有效正/逆势垒分列，不把负差强制改成零。
- DRC是CatMAP官方局部能量敏感度，不等同于全局限速步骤判定。极低TOF体系中的有限差分DRC可能受数值病态影响；诊断表提供来源/日志检查及低TOF提醒。尚无导出的状态Jacobian特征值或多初值吸引域证据，故不把这些稳态结果宣称为动力学可达性证明。
- 机理读法：Ti–Fe 的OHHa覆盖约0.910且OHHa DRC约−0.906，指向对后段OHHa能量/积累敏感；Ti的TS1 DRC约+0.999，是本组中最明显的显式NEB势垒敏感例。多数异核体系（含Ti–W、Ti–Ta）EpoxideHHO接近饱和，且HHO DRC约+1、EpoxideHHO DRC约−1，显示活性响应主要受吸附态能量平衡和中间体释放/转化共同影响；这不是“最高覆盖=唯一限速步”的证明。Ti–Hf、Ti–Ta的TOF极低，DRC仅作警示性局部敏感度，不作强因果结论。

## 13个位点汇总

| 位点 | CatMAP状态 | TOF (s⁻¹) | log10(TOF) | 主导稳定吸附态 | 覆盖度 | 最大|净通量|相对离散 |
|---|---:|---:|---:|---|---:|---:|
| Ti–Fe | CONVERGED | 2.9284154e-05 | -4.53337 | OHHa_s | 0.90974011 | 5.518e-37 |
| Ti–Mn | CONVERGED | 1.5077693e-12 | -11.82167 | EpoxideHHO_s | 1 | 5.834e-26 |
| Ti–Hf | CONVERGED | 6.1710924e-28 | -27.20964 | EpoxideHHO_s | 1 | 2.845e-11 |
| Ti–Re | CONVERGED | 8.3558042e-07 | -6.07801 | EpoxideHHO_s | 1 | 5.005e-63 |
| Ti–Mo | CONVERGED | 2.7918682e-10 | -9.55411 | EpoxideHHO_s | 1 | 1.064e-35 |
| Ti–V | CONVERGED | 7.4403895e-16 | -15.12840 | EpoxideHHO_s | 1 | 2.101e-43 |
| Ti–Co | CONVERGED | 3.7084242 | 0.56919 | EpoxideHHO_s | 0.99962415 | 3.017e-51 |
| Ti–Ti | CONVERGED | 1.0049598e-07 | -6.99785 | EpoxideHHO_s | 0.9961701 | 2.031e-32 |
| Ti–W | CONVERGED | 1.2857614e-15 | -14.89084 | EpoxideHHO_s | 1 | 2.150e-21 |
| Ti–Ta | CONVERGED | 7.2831397e-24 | -23.13768 | EpoxideHHO_s | 1 | 4.531e-20 |
| Ti–Cr | CONVERGED | 1.1506603e-14 | -13.93905 | EpoxideHHO_s | 1 | 2.474e-47 |
| Ti | CONVERGED | 7.6917241 | 0.88602 | H2O_s | 0.0011491636 | 5.874e-67 |
| Ti–Ni | CONVERGED | 0.023417158 | -1.63047 | EpoxideHHO_s | 0.99999977 | 6.248e-39 |

## 各位点最大绝对DRC参数（仅作局部敏感度摘要）

| 位点 | 参数 | DRC | TOF (s^-1) |
|---|---|---:|---:|
| Ti–Fe | OHHa_s | -0.90619 | 2.92842e-05 |
| Ti–Mn | EpoxideHHO_s | -1 | 1.50777e-12 |
| Ti–Hf | EpoxideHHO_s | -1.00251 | 6.17109e-28 |
| Ti–Re | HHO_s | 1.04167 | 8.3558e-07 |
| Ti–Mo | HHO_s | 1.03041 | 2.79187e-10 |
| Ti–V | HHO_s | 1.01637 | 7.44039e-16 |
| Ti–Co | HHO_s | 1.02779 | 3.70842 |
| Ti–Ti | HHO_s | 1.02233 | 1.00496e-07 |
| Ti–W | HHO_s | 1.01486 | 1.28576e-15 |
| Ti–Ta | HHO_s | 1.04595 | 7.28314e-24 |
| Ti–Cr | EpoxideHHO_s | -1 | 1.15066e-14 |
| Ti | TS1_s | 0.998818 | 7.69172 |
| Ti–Ni | HHO_s | 1.02411 | 0.0234172 |


## DRC解释边界

官方D(R)C保留原始符号：正值表示提高该能量在局部会降低产物速率；负值表示相反的局部响应。不同表面态之间存在耦合，单项DRC不等于唯一决速步。对于接近数值下限的TOF，不宜把很大的、微小或符号敏感的有限差分值作强机理断言。具体逐态数据见 `rate_control_matrix.csv` 和 `rate_control_diagnostics.csv`。

## 可复现性与图表

- 复现分析：运行 `scripts/analyze_mechanism_outputs.py`；该脚本只读取CSV、组表和绘图，不调用稳态求解器。
- 补充官方DRC：`scripts/run_missing_official_drc.py`；运行副本与日志保存在本分析目录的 `data/official_drc_runs_supplemental/`。
- H盘目录是归档副本；脚本仍记录/引用本次C盘工作区、CatMAP安装路径及既有网格缓存路径。没有做路径迁移或从H盘独立重放，因此归档副本用于审计保存，不宣称可脱离当前主机目录独立运行。
- PNG/PDF图位于 `../figures/`；Origin矩阵/长表位于 `../origin_exports/`。
- 方向速率由基线CatMAP接口以双精度导出；净通量、覆盖度和TOF保留高精度十进制文本。除格式精度差异外，本数据中第1/6步还出现负方向速率及显著通量闭合异常，根因与官方方向速率路径对近零空位的双精度计算有关；这些方向值不作物理单向通量解释。高精度稳态净通量、覆盖度及TOF的验证独立保留。
