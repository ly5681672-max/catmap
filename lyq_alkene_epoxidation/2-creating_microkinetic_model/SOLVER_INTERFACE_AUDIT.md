# CatMAP 求解器接口与 333.15 K 恢复审计

审计日期：2026-10-11

## 结论

本项目已经切换为只在 333.15 K、有效压力/浓度参数 1.0 下运行的多初值恢复流程。没有继续使用 550 K 或 750 K 的数值辅助路径，也没有修改 DFT 能量、八步机理或反应条件。

9 个原先未通过独立验证的表面均完成了 8 组初值尝试，但没有一个表面满足“至少两个独立初值收敛到同一 TOF 根”的正式提升条件。因此本轮没有向正式基线追加新的催化剂结果。

## CatMAP 接口核对

- 实际使用环境：`G:\anaconda\python.exe` 下的本地 CatMAP；本地包报告 `__version__ = 0.3.1`。
- 官方代码对照：`catmap_official_latest_20261010`，提交 `ed04f91`；其 `numbers_solver.py` 与本地求解器接口一致。
- 实际模型类型：`MinResidMapper` + `SteadyStateSolver`。
- `use_numbers_solver = True`，`numbers_type = None`；CatMAP 默认选择 squared-numbers 转换器。
- 初值向量顺序由 CatMAP 实际返回的 adsorbate 顺序确定，并附加一个 vacancy 项；本项目用
  `x_i = sqrt(theta_i)` 写入 `numbers_map`，再由 CatMAP 执行
  `theta_i = x_i^2 / sum(x_j^2)`。
- 每个初值都记录了 `coverage_vector`、`numbers_vector`、向量长度、位点数、覆盖度和以及往返变换误差，记录文件为 `analysis_ooh/recovery_single_temperature/seed_definitions.csv`。

官方参考：

- [CatMAP numbers solver](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/numbers_solver.py)
- [CatMAP refining a microkinetic model](https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html)

## 反应计量核对

CatMAP 实际解析出 8 个 elementary reactions。多箭头反应使用首态和末态构造净化学计量，而不是把过渡态占位物种当作稳定中间体。解析结果为：

1. `s + H2O2_g -> H2O2_s`
2. `H2O2_s + s -> Ha_s + OOH_s`
3. `s + C6H12_g -> C6H12_s`
4. `OOH_s + C6H12_s + s -> C6H12O_s + O_s + Hb_s`
5. `C6H12O_s -> C6H12O_g + s`
6. `Hb_s + O_s -> OH_s + s`
7. `Ha_s + OH_s -> s + H2O_s`
8. `H2O_s -> H2O_g + s`

`run_ooh_analysis.py` 现已根据 `model.elementary_rxns` 动态构造气相和表面化学计量矩阵，并把以下两项加入独立质量门控：

- 表面物种净生成速率相对于最大步速率的残差；
- CatMAP 报告的气相 turnover frequency 与解析计量矩阵重构值的一致性。

对已完成批处理中的所有 `quality_pass=True` 记录进行回放，气相重构相对误差为 0；表面计量相对误差均低于 `1e-4` 门槛，少数早期精度记录的最大值约为 `2.1e-5`。

## 333.15 K 多初值结果

这里的“有效根”仅表示该初值通过了单根的稳态、符号、覆盖度、跨精度和计量门控；它不表示该根已经是唯一物理稳态。

| 表面 | 有效初值数 | 结果分类 | 代表性 `log10(TOF)` |
|---|---:|---|---:|
| Ti–Mn | 1 | 有效根不足 | -60.4038 |
| Ti–Hf | 1 | 有效根不足 | -120.1156 |
| Ti–Re | 3 | 多个有效稳态分支 | -32.7298、-56.7062 |
| Ti–Mo | 3 | 多个有效稳态分支 | -36.6982、-57.1890 |
| Ti–V | 1 | 有效根不足 | -84.5395 |
| Ti–Zr | 0 | 没有有效根 | — |
| Ti–Ti | 3 | 多个有效稳态分支 | -38.2933、-40.0385、-127.0629 |
| Ti–Ta | 3 | 多个有效稳态分支 | -66.1328、-83.0130 |
| Ti–Cr | 1 | 有效根不足 | -61.6575 |

因此，Ti–Ti 的问题已经从“求解器不收敛”具体化为：在同一 333.15 K 条件下，不同初值会进入不同的数值稳态分支。Ti–Re、Ti–Mo 和 Ti–Ta 也表现出同样现象。对于这些表面，不能按最高 TOF、最低 TOF 或某一个初值结果自动选根。

## 结果文件与基线规则

- 批处理明细：`analysis_ooh/recovery_single_temperature/seed_trials.csv`
- 每个表面汇总：`analysis_ooh/recovery_single_temperature/recovery_results.csv`
- 初值定义与原生向量审计：`analysis_ooh/recovery_single_temperature/seed_definitions.csv`
- 旧的 `recovery_v2` 高温恢复目录保留为历史证据，没有被用于本次新结论。
- 由于本轮没有跨初值一致的表面，`OOH_TOF_results.csv` 没有被追加这些未验证结果，OOH 描述符拟合也没有使用它们。

当前科学结论应写成“333.15 K 下存在初值依赖/多稳态数值分支，剩余表面尚未获得唯一可验证 TOF”，而不是写成某个催化剂已经具有更高活性。尤其不能把 OOH 覆盖度与活性大小规律直接等同。
