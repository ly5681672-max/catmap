# CatMAP 稳态恢复V2：合格种子与多初值验证

## 研究输入保持不变

不修改 `energies_dft_neb.txt`、`DFT_TS_relative.csv`、`alkene_epoxidation.mkm` 中的DFT能量和八步反应。正式TOF计算温度仍为333.15 K，有效储库活度也不变。

## 为什么修正上一版温度路径

CatMAP `MinResidMapper.get_coverage_map` 先检查历史 `numbers_map/coverage_map` 点是否属于新描述符网格。之前从高温缓存复制到仅含333.15 K的单点网格时，高温历史点可能被忽略。V2让 **333.15 K和已验证的高温种子同时位于目标网格**，并用 `run(recalculate=True)` 重新求解。

官方来源：[Refining a Microkinetic Model](https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html)、[MinResidMapper源码](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/mappers/min_resid_mapper.py)、[ReactionModel源码](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/model.py)。逐点质量检查和跨初值TOF一致性是**本项目加严的质量控制**，不是CatMAP默认参数。

## V2 求解方案

1. 对当前未通过验证的9种催化剂建立550/750 K→333.15 K映射（各17点）。
2. 对每一个温度点，用CatMAP已经算出的 `coverage_map`、`rate_map`、`turnover_frequency_map` 核查同一温度下八步净通量、物料守恒和稳态残差。不把内部 `numbers_map` 数值直接当成覆盖度。
3. 只有通过**独立通量及残差校验**的高温点才能打包为新的单点 `coverage_map/numbers_map` 初始种子。
4. 选取最多3个不同温度的合格种子，分别在333.15 K进行260/360/460位精度复算。每级从**相同合格高温初值**重新启动，不传播不合格目标解。
5. 每个初值要求连续两个精度结果均通过且TOF一致；至少两个独立温度种子得出的 `|Δlog10(TOF)| ≤ 0.02` 才可更新最终基准表。存在分支分歧则记录，不擅自选择某一根。
6. 原先已经数值验证的Ti–Fe、Ti–Co、Ti–Ni、Ti–W、Ti不会被未通过的恢复结果覆盖。

这些阈值是数值质量控制，并非对稳态解唯一性或实验活性的证明。

## Windows运行

```powershell
cd H:\catmap\catmap
git -c http.proxy=http://127.0.0.1:6696 -c https.proxy=http://127.0.0.1:6696 pull --ff-only origin main
cd lyq_alkene_epoxidation\2-creating_microkinetic_model
python -m unittest -v test_ooh_numerics test_ooh_recovery
python run_ooh_analysis.py --recover --recover-surfaces titi timn ticr
python run_ooh_analysis.py --recover
```

## 新版输出

- `analysis_ooh/recovery_v2/bridge_diagnostics.csv`：每条温度路径的映射点数及质量合格种子数。
- `analysis_ooh/recovery_v2/bridge_point_quality.csv`：每个温度点的通量和残差验证记录。
- `analysis_ooh/recovery_v2/seed_trials.csv`：每个高温种子独立精度计算的详细记录。
- `analysis_ooh/recovery_v2/recovery_results.csv`：各催化剂跨种子一致性及最终判定。

新恢复模块**不覆盖旧版 `analysis_ooh/recovery/` 结果**。只有额外验证通过的催化剂才会更新 `dft_baseline.csv`、`OOH_TOF_results.csv` 和相应拟合图。GitHub代码提交不等于已在本地完成新的CatMAP真实求解，需运行后查看V2 CSV。
