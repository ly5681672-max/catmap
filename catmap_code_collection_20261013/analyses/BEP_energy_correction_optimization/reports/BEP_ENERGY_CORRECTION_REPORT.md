# Ti–M氧转移NEB1扩展BEP与能量修正分析

## 研究边界与数据
本分析仅使用已归档的13个位点严格CatMAP形成能、历史单吸附态OH*/Ha*/OOH*数值及严格模型OHHa态；未运行DFT、CatMAP稳态/TOF或条件扫描。正式目标为有符号 `Y = G(TS1*) − G(Rpair*)`，传统BEP自变量为 `ΔG₂ = G(EpoxideHHO*) − G(Rpair*)`。负的TS1−Rpair值保留。Ti–Zr不纳入。历史OH*/Ha*/OOH*参考的slab闭合尚未完成，相关扩展系数均属探索性统计关系，不解释为已验证因果效应。

样本按surface_name一对一匹配；13个正式位点为：Ti-Fe, Ti-Mn, Ti-Hf, Ti-Re, Ti-Mo, Ti-V, Ti-Co, Ti-Ti, Ti-W, Ti-Ta, Ti-Cr, Ti, Ti-Ni。11位点对照仅排除Ti-W、Ti-Ta，未替代正式13位点结果。

## 模型与复现指标
OLS系数按eV输入与eV输出；LOOCV每次仅用训练折重新估计系数。11-site refit与13-site系数用于11点是两类不同评估。

|模型|13-site R²|RMSE (eV)|MAE (eV)|LOOCV RMSE (eV)|11-site refit R²|11-site refit RMSE|13-site系数用于11点 R²|
|---|---:|---:|---:|---:|---:|---:|---:|
|Traditional BEP|0.2054|0.4177|0.3060|0.7134|0.2935|0.1293|-1.4240|
|BEP + OH + Ha + OOH|0.6306|0.2848|0.2226|0.7136|0.7270|0.0804|-1.4742|
|BEP + OH + Ha + OHHa|0.6757|0.2669|0.2149|0.6710|0.7042|0.0837|-1.2047|
|OH-modulated slope|0.7564|0.2313|0.1743|0.4197|0.3103|0.1278|-0.7294|

### 13位点拟合公式

- **Traditional BEP**：`Y = (+0.367197)·deltaG2_eV + (+0.907242 eV)`
- **BEP + OH + Ha + OOH**：`Y = (-0.056897)·deltaG2_eV + (+0.311555)·OH_star_legacy_eV + (-0.22131)·Ha_star_legacy_eV + (-0.00535239)·OOH_star_legacy_eV + (+0.949674 eV)`
- **BEP + OH + Ha + OHHa**：`Y = (-0.199181)·deltaG2_eV + (+0.323272)·OH_star_legacy_eV + (-0.317337)·Ha_star_legacy_eV + (+0.256458)·OHHa_strict_eV + (+2.4956 eV)`
- **OH-modulated slope**：`Y = (+0.855092)·deltaG2_eV + (-0.432901)·deltaG2_x_OH_centered + (-0.75469)·OH_centered_eV + (+2.06922 eV)`

## 结果要点
13位点样本内最高R²为OH变斜率模型（0.7564），其LOOCV RMSE=0.420 eV，较传统BEP低；但11位点重拟合R²仅0.310、LOOCV RMSE=0.614 eV，显示模型对样本构成敏感。方案B在11位点的LOOCV RMSE最低（0.142 eV），但将13点系数直接应用于11点时R²为-1.205，说明其跨子集迁移也没有得到支持。
嵌套消融显示OH与Ha加入后明显改善样本内拟合，A中的OOH追加贡献近零；B中OHHa有一定增量。完整模型中ΔG₂系数反号且删去ΔG₂后的R²损失较小，因此改善主要不是传统正斜率BEP项得到强化，而是能量修正项/交互结构对当前样本的解释增加。历史OH/Ha/OOH的参考不确定性、变量共线性和小样本仍限制物理解读。

## 变斜率定义
OH中心值为 `-4.2305517 eV`。模型为 `Y = α₀ΔG₂ + α₁ΔG₂·Z + γZ + β`，其中 `Z = E_OH − mean(E_OH)`；13位点 `α₀=0.855092`，`α₁=-0.432901`，有效斜率 `α_eff=α₀+α₁Z` 的样本范围为 `-0.00466345` 至 `1.83077`。此为样本内交互模型，不等同已验证的机制斜率定律。

## 贡献与稳健性诊断
完整项与逐项删除、ΔG₂删除、传统BEP斜率固定后的结果见 `BEP_CONTRIBUTION_DIAGNOSTICS.csv`；变量相关、VIF与标准化设计矩阵条件数见 `BEP_COLLINEARITY_DIAGNOSTICS.csv`；逐一留出一个催化剂重拟合的系数变化见 `BEP_LEAVE_ONE_COEFFICIENT_SENSITIVITY.csv`。这些量用于诊断，不将多变量OLS系数当作独立因果贡献。
按删除项后的R²差衡量，方案A去掉ΔG₂仅下降-0.0029；方案B去掉ΔG₂下降-0.0254。相反，把ΔG₂斜率固定为传统BEP的13点斜率后，A/B的R²分别为0.4686/0.4699。A中OOH项的条件增量R²仅约0.00003；因此不能把A的提升归因于OOH单项。

最大绝对样本内残差（按模型）：model_A_OH_Ha_OOH: Ti-Ta (-0.6004 eV)；model_B_OH_Ha_OHHa: Ti-Ta (-0.5236 eV)；traditional_BEP: Ti-Ta (-1.1813 eV)；variable_slope_OH: Ti-Cr (+0.5155 eV)

Ti–W、Ti–Ta、Ti–Co、Ti–Fe的单点残差在 `BEP_CATALYST_RESIDUALS.csv` 保留。Ti–W/Ti–Ta的有符号NEB1差值为负值，不做截断。

|催化剂|传统BEP残差|方案A残差|方案B残差|OH变斜率残差|
|---|---:|---:|---:|---:|
|Ti-Fe|+0.0523|-0.1717|+0.0839|-0.1323|
|Ti-Co|-0.2629|-0.1920|-0.1806|-0.0151|
|Ti-W|-0.4914|-0.2233|-0.2786|+0.1956|
|Ti-Ta|-1.1813|-0.6004|-0.5236|-0.4548|
逐一移除一个位点后重拟合，并比较该系数集对全部13位点的预测变化，影响最大的留出位点分别为：Traditional BEP: tico（RMS预测变化 0.507 eV）；BEP + OH + Ha + OOH: tico（RMS预测变化 0.533 eV）；BEP + OH + Ha + OHHa: tico（RMS预测变化 0.488 eV）；OH-modulated slope: tico（RMS预测变化 0.221 eV）。这是小样本影响诊断，不等价于删除这些点。

### 逐项消融（13-site）

|模型序列|新增项|R²|RMSE (eV)|
|---|---|---:|---:|
|A_OH_Ha_OOH|deltaG2_eV|0.2054|0.4177|
|A_OH_Ha_OOH|OH_star_legacy_eV|0.4934|0.3335|
|A_OH_Ha_OOH|Ha_star_legacy_eV|0.6305|0.2848|
|A_OH_Ha_OOH|OOH_star_legacy_eV|0.6306|0.2848|
|B_OH_Ha_OHHa|deltaG2_eV|0.2054|0.4177|
|B_OH_Ha_OHHa|OH_star_legacy_eV|0.4934|0.3335|
|B_OH_Ha_OHHa|Ha_star_legacy_eV|0.6305|0.2848|
|B_OH_Ha_OHHa|OHHa_strict_eV|0.6757|0.2669|

## TS1重建
每种模型以 `G(TS1)_pred = G(Rpair)_actual + Y_pred` 重建。`BEP_TS1_RECONSTRUCTION.csv` 同时包含样本内与LOOCV预测及残差。这里仅诊断能量重建，不重算TOF。

## 图件与可复现文件
7类图均输出PNG/PDF：传统BEP、方案A预测一致性、方案B预测一致性、四模型残差、OH调制斜率机制图、13/11拟合性能对照、变量消融。图中保留13个位点名称，Ti–W/Ti–Ta使用不同标记。Origin所需原始点、拟合值、残差、系数和曲线数据位于 `origin_exports/`。

复现：在分析包根目录运行 `python scripts/run_bep_energy_correction.py`。脚本只读取 `data/` 快照并写入本包，不调用CatMAP运行接口。

## 解释限制
训练R²提高不保证外推；13个样本而多变量项数较多，LOOCV及逐位点系数变化需一并考虑。历史单吸附态数据与严格共吸附模型能量的slab/热修正参考尚未闭合，因此方案A中的OH*/Ha*/OOH*及交叉比较具有参考不确定性。方案B使用严格能量表OHHa态，但仍是跨反应网络能量项的统计耦合，不能据此断言OHHa直接控制NEB1。

## 参考与软件
数值拟合使用本机NumPy OLS（`numpy.linalg.lstsq`）、逐点留一重拟合及Matplotlib；CatMAP官方结果仅作为只读数据源，本分析未调用 `ReactionModel.run()`。
图件设计采用原始数据与派生分析分离、图例/符号双重编码、保留全部点与有符号值等原则；参考：Kassis T, Agarwal V, He Y, Patel D, Brueckner AM. Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065v2 (2026). https://arxiv.org/abs/2609.00065
