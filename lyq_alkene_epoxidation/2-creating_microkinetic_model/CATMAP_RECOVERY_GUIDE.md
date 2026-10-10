# CatMAP官方求解路径恢复（保持原始DFT不变）

## 本轮目的

目标是恢复当前14种催化剂中**尚未通过**跨精度稳态验证的9种，不修改：
- \`energies_dft_neb.txt\` 等DFT能量；
- \`DFT_TS_relative.csv\`；
- 八步反应 \`rxn_expressions\`；
- 333.15 K的正式反应条件以及原有的气相有效储库活度。

已通过的Ti–Fe、Ti–Co、Ti–Ni、Ti–W和Ti不会被重新覆盖。

## 官方机制与对应源码

1. [CatMAP: Refining a Microkinetic Model](https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html)明确说明可以使用**之前求解的覆盖度作为下一次模型的初始猜测**。这种方式可以提高刚性模型的收敛概率。
2. [CatMAP MinResidMapper源代码](https://catmap.readthedocs.io/en/latest/_modules/catmap/mappers/min_resid_mapper.html)说明映射器从邻近描述符点寻找初始解，在困难路径中按\`max_bisections\`对描述符距离二分，且\`use_numbers_solver=True\`时可读取已保存的\`numbers_map\`。
3. [CatMAP ReactionModel.run源码](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/model.py)支持\`run(recalculate=True)\`：**用历史数据作为初始猜测，但是重新求解**，避免把之前的结果误当成新收敛结果。
4. 使用\`ThermodynamicScaler\`的\`temperature\`、\`pressure\`描述符，临时建立高温→333.15 K的温度路径。温度路径属于**本项目增加的数值homotopy策略**，不是CatMAP官方规定的标准温度；正式TOF始终在333.15 K重新计算。

## 执行机制

\`\`\`text
保持原DFT能量 + 原八步反应
          │
          ▼
任选550 K / 750 K作为数值路径高温端
          │
          ▼
MinResidMapper沿温度网格逐步逼近333.15 K
  （邻近覆盖度 + max_bisections=5）
          │
          ▼
保存CatMAP原生coverage_map和numbers_map
          │
          ▼
在333.15 K目标点用缓存作为初始猜测
          │
          ▼
260/360/460位逐级复算，严格独立校验
          │
          ▼
连续两级验证通过，才更新dft_baseline和OOH–TOF
  （否则记录错误，不修改原有通过结果）
\`\`\`

注意：即使温度路径能找到某个数值根，也不能保证它是物理上唯一稳定的稳态解；仍以基元通量、残差、覆盖度范围、气相守恒及跨精度TOF一致性为硬性校验。

## 使用方式（Windows，需本地安装CatMAP）

\`\`\`powershell
cd H:\catmap\catmap
git -c http.proxy=http://127.0.0.1:6696 -c https.proxy=http://127.0.0.1:6696 pull --ff-only origin main

cd H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model
python -m unittest -v test_ooh_numerics test_ooh_recovery

# 优先诊断3种第四级失败的催化剂
python run_ooh_analysis.py --recover --recover-surfaces titi timn ticr

# 对目前仍未通过验证的其余体系执行恢复
python run_ooh_analysis.py --recover

# 如有需要可增加温度路径选择，但并不意味着修改目标反应温度
python run_ooh_analysis.py --recover --bridge-temperatures 450 600 800
\`\`\`

## 输出结果

- \`analysis_ooh/recovery/bridge_diagnostics.csv\`：每个表面的桥接温度、已成功求解的温度点数、是否覆盖333.15 K；
- \`analysis_ooh/recovery/target_refinement.csv\`：每级目标温度333.15 K的八步净反应通量、残差与TOF；
- \`analysis_ooh/recovery/recovery_results.csv\`：每个待修复表面的最终状态；
- \`analysis_ooh/recovery/baseline_before_recovery.csv\`：首次成功恢复时自动保留的原基准表；
- \`analysis_ooh/dft_baseline.csv\`：**仅新增通过跨精度验证的催化剂行**；
- \`analysis_ooh/OOH_TOF_results.csv\`及\`OOH_logTOF_linear.png/.pdf\`：若有新通过项则基于更新后的合格数据重建。

最重要的规则：**未找到有效解不等于催化剂没有活性；只表示本轮CatMAP数值恢复仍不充分。** 原始DFT可信性和数值求解的通过与否是两个不同维度。

## 测试与验证边界

\`test_ooh_recovery.py\`只检查恢复元数据、目标温度识别和参数输入，不会假装本地CatMAP已经完成真实温度路径求解。温度映射、缓存加载、真实求解与数值稳定性，需要按上述命令在用户的CatMAP环境运行后检查CSV和日志。
