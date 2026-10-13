# O*/Ha*/Hb* 历史形成能—严格 CatMAP TOF 经验图

## 数据及图的性质

- 位点：Ti-Fe, Ti-Mn, Ti-Hf, Ti-Re, Ti-Mo, Ti-V, Ti-Co, Ti-Ti, Ti-W, Ti-Ta, Ti-Cr, Ti, Ti-Ni（13个；Ti–Zr未收敛而排除，Ti–Nb不属于本轮严格基准集合）。按 `surface_name` 一对一连接；没有按行号拼接。
- y 数据：各位点 `official_strict_tolerance/<surface>/gas_turnover_frequencies.csv` 中 `C6H12O_g`，并核对 `convergence.csv` 为 CatMAP 0.3.1、333.15 K、`CONVERGED`。未重算 CatMAP。
- x 数据：原有 `DFT_formation_catmap.txt` 中 O、Ha、Hb `formation_energy`，其转换脚本使用 `G_ads − G_slab − μ_atom`，其中旧单物种数据分别来自 `catmap数据集` 的 Gibbs 缓存列，slab 为第51行，原子参考为 μH=-1.11 eV、μO=-4.37 eV；Ha 与 Hb 独立保留。复算 39 个值全部与旧表一致。
- 每张图上的圆点/菱形及名称是实测严格 CatMAP 位点数据；彩色曲线和二维色面是对这些点做的经验 OLS 拟合，不是 CatMAP 网格计算、不是新的 TOF。二维拟合限在 13 点描述符凸包内。
- 一维图同时显示线性、二次、连续非对称分段线性拟合；二维图显示 O+Ha 或 O+Hb 六参数二次响应面。峰或最大值均由对应拟合函数在样本范围/凸包中计算；TOF=10^{log10(TOF)}。

## 拟合结论与峰值

| 描述符 | 样本内二次 R² | 线性/二次/分段 LOOCV RMSE（log10 TOF） | 峰判断 |
|---|---:|---:|---|
| O* | 0.419 | 8.63 / 7.55 / 8.61 | 无内部火山峰；分段斜率为 -1.18, 8.82，呈谷形转向。二次模型最大值在边界。 |
| Ha* | 0.246 | 11.35 / 24.97 / 18.33 | 无内部火山峰；分段斜率为 -7.25, 6.90，呈谷形转向。二次模型最大值在边界。 |
| Hb* | 0.795 | 5.29 / 4.79 / 4.68 | 二次模型有样本范围内峰：Hb=-4.114 eV，预测 log10TOF=-0.731，TOF=0.186；分段模型两侧斜率均为负，不支持火山。 |

二维 O+Ha：样本内 R²=0.679，LOOCV RMSE=16.95；最大值为凸包边界最大值，O=1.980 eV、Ha=-2.121 eV，预测 log10TOF=2.025、TOF=106，邻近真实点：Ti-Co, Ti-Fe, Ti-Ti。

二维 O+Hb：样本内 R²=0.896，LOOCV RMSE=4.82；最大值为凸包边界最大值，O=1.922 eV、Hb=-4.475 eV，预测 log10TOF=1.132、TOF=13.6，邻近真实点：Ti, Ti-Ni, Ti-Fe。这不是内部二维最优活性区。

## 能量参考与温度证据

旧表与当前严格网络不是已验证的统一能量参考。旧转换器用 `catmap数据集` 每个表面第51行作为 clean slab；严格网络用 `吉布斯自由能汇总` 每个表面第4行作为 clean slab。当前工作簿两处 slab 数值差异（current−old）≥0.05 eV 或 row41−row4 差异≥0.05 eV 的位点：

- Ti-Fe: old slab -870.962186 eV (G51); current row-4 slab -869.897259 eV (F4), Δ(current−old)=+1.064927 eV; row41−row4=-1.063296 eV.
- Ti-Mo: old slab -873.868832 eV (AS51); current row-4 slab -873.920517 eV (AR4), Δ(current−old)=-0.051685 eV; row41−row4=+0.059830 eV.
- Ti-Co: old slab -868.516201 eV (BN51); current row-4 slab -867.638499 eV (BM4), Δ(current−old)=+0.877702 eV; row41−row4=-0.000512 eV.
- Ti-Ti: old slab -873.864081 eV (BU51); current row-4 slab -873.864076 eV (BT4), Δ(current−old)=+0.000005 eV; row41−row4=-0.766641 eV.
- Ti-Ta: old slab -877.703381 eV (CI51); current row-4 slab -877.574536 eV (CI4), Δ(current−old)=+0.128845 eV; row41−row4=-0.128640 eV.
- Ti-Ni: old slab -863.139429 eV (DE51); current row-4 slab -863.254380 eV (DD4), Δ(current−old)=-0.114951 eV; row41−row4=+0.117033 eV.

旧单物种 O/Ha/Hb 的 `T` 列（catmap数据集对应各源行）对 39/39 记录为空；虽然能追溯其 `E0`、`ΔG` 修正列及缓存公式，但没有证据把这些修正归属到 333.15 K。严格汇总表第4行的温度字段可记录 clean-slab 项的温度，但不能补证旧 O/Ha/Hb 单物种热修正。旧与新 slab 差异不得自动作为校正量，因为结构等价性和热力学基准尚未验证。故**不生成所谓统一参考 O/Ha/Hb 表**，本轮图只作为历史参考下的探索性经验关系。逐单元格、公式、缓存值和 slab 对照见 `workbook_source_audit.csv`、`reference_slab_comparison.csv`；输入哈希见 `fit_metrics.json`。

## 文件

- `O_TOF_empirical.png/.pdf`
- `Ha_TOF_empirical.png/.pdf`
- `Hb_TOF_empirical.png/.pdf`
- `O_Ha_TOF_empirical_map.png/.pdf`
- `O_Hb_TOF_empirical_map.png/.pdf`
- `matched_data.csv`：逐 surface_name 匹配的13条严格 TOF 和旧形成能。
- `fit_metrics.json`：拟合系数/指标、峰值及源文件 SHA256。
- `workbook_source_audit.csv`：39个物种源单元格公式、缓存值、ΔG/E0、温度字段与重算核验。
- `reference_slab_comparison.csv`：旧 slab、新 row4 slab 与 row41 参照差异。
- `analyze_legacy_descriptors.py`：只读数据的可复现绘图脚本。
