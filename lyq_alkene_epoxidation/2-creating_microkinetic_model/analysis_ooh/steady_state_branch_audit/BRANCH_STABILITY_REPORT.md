# CatMAP 333.15 K 稳态分支与 TOF 审计报告

日期：2026-10-11。研究体系：14 种 Ti/Ti–M 催化剂的 1-己烯/H₂O₂ 环氧化。所有本轮数值求解固定在 `[333.15 K, 1.0]`；原始 DFT、TS 映射和八步反应式的 SHA256 与 `input_manifest.json` 一致。既有正式 TOF 表没有更新。

## 结论

本模型的八步网络除位点守恒外，还精确保留两个独立的表面库存：

\[
I_1=\theta_{\mathrm{Ha}}-\theta_{\mathrm{OOH}}-\theta_{\mathrm{OH}}-\theta_{\mathrm O},
\qquad I_2=\theta_{\mathrm{Hb}}-\theta_{\mathrm O}.
\]

这不是 CatMAP 求解器故障，而是从 CatMAP 实际解析的八步反应矩阵算出的精确零空间：9 个吸附物种的化学计量矩阵秩为 7，零空间维数为 2。`stoichiometric_invariants.json` 保存了物种顺序、解析后的反应和核验结果。若真实动力学从洁净空表面出发，始终有 `I1=I2=0`；如催化剂经预吸附/预处理，则须由实验或明确的初始覆盖度确定这两个值。

79 个通过跨精度、独立速率重算、残差和计量检查的数值根，按完整覆盖度与八步净速率聚为 43 组不同的**采样稳态解**。43 个代表根在各自固定 `(I1,I2)` 的 7 维动力学子空间内，Jacobian 特征值实部均为负且导数交叉检查通过；但在完整物理覆盖度单纯形内均有两个守恒零模，故完整空间标签是 `near_neutral_or_undetermined`。它们不能称作 43 个彼此竞争的、全空间渐近稳定吸引子。与空表面零库存相容的根仅有 11 个：Ti–W 10 个，Ti–Ti 1 个。

目前没有证据支持将 14 种催化剂都赋予唯一可靠 TOF，也不应以任何非零库存数值根的最大/最小 TOF 进行活性排序或 OOH*–TOF 拟合。Ti–W 的条件性候选 TOF 为 `1.225349870552975109149936162064160132845e-137`（`log10 TOF=-136.91173989076393`，沿用项目 TOF 单位）；这仅在零库存初始类及现有模型成立时有支持，不自动更新正式结果。Ti–Ti 的零库存候选为 `9.151047062590976e-41`，仅一个合格根，仍需独立复现。

## 14 种催化剂统一判定

“簇”表示有限初值扫描到的不同完整覆盖度数值解；它们可能属于不同守恒库存，不能直接解释为同一实验条件下的多重吸引子。`0 类根`列只统计 `|I1|,|I2| ≤ 1e-10` 的候选；该容差是数值分组标准，不是修改物理模型。

| 催化剂 | 合格根 | 采样簇 | 0 类根 | 零库存条件下结论 |
|---|---:|---:|---:|---|
| 单 Ti | 12 | 5 | 0 | 未找到合格根 |
| Ti–Co | 8 | 3 | 0 | 未找到合格根 |
| Ti–Cr | 4 | 3 | 0 | 未找到合格根 |
| Ti–Fe | 10 | 5 | 0 | 未找到合格根 |
| Ti–Hf | 2 | 1 | 0 | 未找到合格根 |
| Ti–Mn | 3 | 2 | 0 | 未找到合格根 |
| Ti–Mo | 3 | 3 | 0 | 未找到合格根 |
| Ti–Ni | 13 | 7 | 0 | 未找到合格根 |
| Ti–Re | 3 | 3 | 0 | 未找到合格根 |
| Ti–Ta | 3 | 3 | 0 | 未找到合格根 |
| Ti–Ti | 3 | 3 | 1 | 单根，证据不足 |
| Ti–V | 3 | 2 | 0 | 未找到合格根 |
| Ti–W | 10 | 1 | 10 | 单簇，零库存条件下有支持 |
| Ti–Zr | 2 | 2 | 0 | 未找到合格根 |

“未找到”是有限初值求解结果，不是数学上的无根证明，更不能写成 TOF 为零。Ti–W 的 10 个根由不同求根初值获得，复现了完整覆盖度和速率；求根初值不是动力学初始条件，因此尚未构成轨迹吸引域证明。`surface_branch_summary.csv` 给出机器可读的分类和正式 TOF 动作，`branch_inventory_summary.csv` 列出每个采样簇的覆盖度主导物种、TOF、守恒库存及两种稳定性标签。

## 重点体系

- Ti–Re：Hb* 主导与 O* 主导根的 `log10 TOF` 均约 `-32.730`，但覆盖度向量显著不同，且 `(I1,I2)` 分别近似 `(0,1)` 与 `(-1,-1)`。独立重算与重放没有发现速率缓存错配；相同 TOF 不能合并为同一根，也不能据此证明同一守恒类中的双稳态。
- Ti–Mo：Hb* 与 O* 主导根均约 `-36.698`，对应上述不同库存；另有 Ha* 主导根约 `-57.189`。三者都是不同采样解。
- Ti–Ta：Hb* 与 O* 主导根均约 `-66.133`，库存分别近似 `(0,1)` 与 `(-1,-1)`；Ha* 主导根约 `-83.013`。
- Ti–Ti：三个候选根的 `log10 TOF` 约 `-127.063`、`-38.293`、`-40.039`。仅最后一个与零库存相容，且只有一个合格求根初值；不能将另外两个直接用于空表面条件下的 TOF 选择。
- Ti–Zr：扩展初值后新增两个跨精度合格根，故原先“完全无有效根”已不适用；但两者 `I1≈-1`，均不对应空表面零库存。空表面初值的三个精度阶段虽返回 `solved`，表面相对稳态残差约为 `2`、解析表面平衡相对误差为 `1`，因此被质量门控拒绝。零库存根未得到的根本原因仍未确定。

## 数值方法与证据边界

1. 原始 16 个跨精度候选先在独立目录重新实例化 CatMAP，清理速率缓存，并重新计算完整覆盖度、numbers、八步净速率、气相 TOF 与动态计量残差。扩展批次的 63 个候选保留了相同的完整状态和清缓存重算检查。`combined_validated_roots.csv` 保存高精度字符串和源目录；`combined_root_clusters.csv` 以完整向量聚类，而非仅凭 TOF 或主导物种。
2. 对每个不同采样簇检查 CatMAP `ideal_steady_state_function` 与 `ideal_steady_state_jacobian` 的物理覆盖度导数，并以约束单纯形内有限差分复核。完整 9 维吸附覆盖度 Jacobian 含两个守恒零模；投影到秩 7 的固定库存子空间后全部 43 个代表根为 `linearly_stable`。`combined_jacobian_eigenvalues.csv` 保留两组特征值、条件数、精度和导数误差。原 16 根的先前 Jacobian 文件作为历史证据保留，最终标签已在合并文件中修正。
3. 本轮没有宣称完成刚性动力学轨迹吸引域验证。严格守恒量已经证明不同库存间的真实轨迹不能相互到达；极低速率体系的有限时间双精度积分也不足以证明趋于稳态。对接近同一库存而 TOF 不同的采样根，是否在**同一精确库存**存在多个吸引子仍待约束动力学求解。
4. CatMAP 官方文档说明 `SteadyStateSolver` 的 tolerance 是表面覆盖度变化率的容许阈值，精度设置用于处理刚性方程；这不等价于证明动力学吸引域或实验初始库存。[CatMAP 官方教程](https://catmap.readthedocs.io/en/latest/tutorials/creating_a_microkinetic_model.html)、[官方稳态求解器源码](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/steady_state_solver.py)。实际接口和结论以本机 CatMAP 0.3.1 的解析反应与输出为准。

## 后续科研决策

先确认实验前催化剂是否可近似看作洁净空表面，或给出预吸附后的 `I1,I2` 初始库存。确定库存后，可在**固定 333.15 K、原八步机理与原 DFT 输入**下增加显式守恒约束的高精度稳态求解，并对同一库存的候选做刚性轨迹或其他吸引域核验。在此之前，正式 TOF 表保持历史状态；OOH*–TOF 拟合不应混入库存未定义或未验证的数值根。若以后需要移除这两个额外守恒量，必须以有物理证据的新反应步骤重新建模，那属于另一项科研决策，本轮未改动网络。

## 文件、测试和可复核性

- 新增 `finalize_branch_audit.py`、`audit_conserved_inventory.py`；修正 `branch_audit.py` 的全覆盖度稳定性标签；`test_branch_audit.py` 增加精确守恒与零模回归测试。
- 审计目录新增 `combined_all_seed_trials.csv`、`combined_validated_roots.csv`、`combined_root_clusters.csv`、`combined_jacobian_eigenvalues.csv`、`conserved_inventory.csv`、`stoichiometric_invariants.json`、`branch_inventory_summary.csv`、`surface_branch_summary.csv`、`numerical_failures.csv`。旧 `recovery_single_temperature/`、正式 TOF 表和原 16 根审计表仍保留。
- `G:\anaconda\python.exe -m unittest discover -p 'test_*.py' -v`：26 项通过。原始 `alkene_epoxidation.mkm`、`energies_dft_neb.txt`、`DFT_TS_relative.csv` 的 SHA256 均与审计开始时一致。

科研证据审查方法参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. [arXiv:2609.00065](https://arxiv.org/abs/2609.00065)。本报告的化学计量推导与数值结果来自本项目 CatMAP 文件，并非该方法论文的研究结果。
