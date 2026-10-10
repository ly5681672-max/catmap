# CatMAP OOH/TOF 数值稳态复算方案（2026-10-10）

## 前提：原始DFT计算数据保持完全不变

按本项目研究设定，\`energies_dft_neb.txt\`、\`DFT_TS_relative.csv\` 和其余DFT原始数据是已确认的输入。**此工作流不会编辑原始DFT数值、改变八步反应网络或以某项相对TS能量为负来自动排除催化剂。**

目的仅是让CatMAP求解器可靠求出在这些已给定能量、反应网络和动力学假设下的稳态通量，并评估OOH*描述符与TOF的关系。

## 1. CatMAP官方依据

- [Creating a Microkinetic Model](https://catmap.readthedocs.io/en/latest/tutorials/creating_a_microkinetic_model.html)：对于刚性反应，使用多精度算法、合理的\`tolerance\`及\`max_rootfinding_iterations\`。
- [Refining a Microkinetic Model — Refining Numerical Accuracy](https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html)：可通过逐步提高\`decimal_precision\`、收紧\`tolerance\`，比较结果对数值参数的独立性；\`threshold\`是**绘图噪声处理**而不是物理TOF界限。
- [CatMAP numbers_solver.py](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/numbers_solver.py)：\`use_numbers_solver=True\`时收敛范数是**表面稳态残差的L2平方和**，因此容差\`1e-120\`约对应原始残差\`1e-60\`。真实TOF若更小，单次绝对收敛并不足以验证相对通量。

## 2. 程序改进

\`run_ooh_analysis.py\` 采用保留原始DFT能量的\`ThermodynamicScaler\`，并对每种催化剂进行独立复算：

| 精度级 | decimal_precision | tolerance | 备注 |
|---|---:|---:|---|
| 1 | 180 | 1e-120 | 初步求解 |
| 2 | 260 | 1e-180 | 收紧求解 |
| 3 | 360 | 1e-260 | 默认最大精度 |
| 4（可选） | 460 | 1e-340 | 对超小通量追加复算 |

注意：所有容差都由\`mpmath.mpf\`精确设置。若\`1e-340\`以原生Python float字面量赋值会下溢到0。

每个精度级保存在各自的\`analysis_ooh/refinement/stage_XX/<surface>/\`目录，避免旧缓存自动充当新结果。

## 3. “有效TOF”的数值定义（本项目额外增加的质量门槛）

按原有8步单一闭合反应循环：对8个净基元速率\`r_1...r_8\`，检查

- 相对通量误差：\`max(|r_i-r_4|) / max(|r_i|) < 1e-4\`；
- 气相计量误差：\`max(|r_C6H12O, -r_C6H12, -r_H2O2, r_H2O - r_C6H12O|) / max(abs(net flux)) < 1e-4\`（以脚本的逐项差值为准）；
- 独立表面稳态残差（numbers_solver取平方根）/最大基元通量 < 1e-4；
- 覆盖度总和、符号、绝对残差合理；
- **连续两级独立求解均通过上述检查，且\`|Δlog10(TOF)| ≤ 0.02\`** 才标记为\`validated_refinement_stable\`。

这些相对阈值及\`0.02\`对数一致性门槛是本项目的**质量控制选择，不是CatMAP官方默认参数**。

**取消固定\`TOF < 1e-40\`判为无效**的旧规则，改为“相对数值残差+连续两次稳定复算”。超小TOF保留\`net_C6H12O_high_precision\`原始高精度字符串和\`log10_net_C6H12O\`，避免转换成float后下溢。

\`TS_requires_check\`字段现在仅保留为原表描述性信息，**不会排除DFT已确认的催化剂**。

## 4. 运行流程（Windows / Python环境需已安装CatMAP）

在\`H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\`中执行：

\`\`\`powershell
python -m unittest -v test_ooh_numerics
python run_ooh_analysis.py --diagnose
python run_ooh_analysis.py --diagnose --coverage-crosscheck
python run_ooh_analysis.py
# 仅在默认3级不能解决且有计算资源时：
python run_ooh_analysis.py --max-stage 4
\`\`\`

\`--diagnose\`专门比较Ti–Fe、Ti–Ti、Ti–W新旧参数，输出\`analysis_ooh/diagnostics/diagnostics_comparison.csv\`；\`--coverage-crosscheck\`额外尝试传统coverage求解器。原日志显示其因Jacobian数值奇异而失败，这不是DFT错误的直接证据。失败将被记录，不能假装输出有效TOF。

14种催化剂的复算输出：

- \`analysis_ooh/refinement_convergence.csv\`：逐级精度和各步通量、独立误差、TOF一致性；
- \`analysis_ooh/dft_baseline.csv\`：每个催化剂最终质量状态；
- \`analysis_ooh/OOH_TOF_results.csv\`：OOH、OH、O形成能、TOF、\`eligible_for_tof_fit\`；
- \`analysis_ooh/OOH_OH_linear.png/.pdf\`：始终保留原始DFT形成能标度分析；
- \`analysis_ooh/OOH_logTOF_linear.png/.pdf\`：仅当至少5种催化剂通过跨精度验证时输出，否则删除本地失效旧图，避免误用。

旧的CSV数据如果没有\`refinement_stable\`及\`log10_net_C6H12O\`字段，\`--fit-only\`会全部判为未验证。

## 5. 当前保持不变的动力学前提

- 八步基元反应与所有原始DFT能量不变；
- \`temperature=333.15 K\`；
- 保留原脚本的\`pressure_mode='concentration'\`及原有有效储库浓度/活度近似，\`C6H12_g=1.0\`、\`H2O2_g=0.185\`、\`H2O_g=0.815\`、\`C6H12O_g=1e-20\`；
- \`gas_thermo_mode='frozen_gas'\`、\`adsorbate_thermo_mode='frozen_adsorbate'\`；
- 原配置没有显式\`prefactor_list\`，沿用CatMAP默认前因子，不擅自引入新的实验动力学参数。

这些是**模型约定**，不代表用计算收敛就能证明液相活度模型、反应机理及前因子与实验完全一致。新的CSV包含条件审计字段，帮助后续敏感性研究。

**重要：GitHub上已有历史数据是旧算法输出；提交新代码不代表重新运行。必须先在本地真实CatMAP环境运行后才判断新TOF和拟合结果。**
