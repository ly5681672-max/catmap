# O*+H*经验二次响应面的全矩形网格导出

## 方法与数据边界
- 输入：本目录 `matched_data.csv`（13 个按 surface_name 匹配的已验证 CatMAP TOF 与历史 O/Ha/Hb 形成能）及 `fit_metrics.json` 中保存的原二次响应面系数。
- 沿用原六参数标准化二次响应面；未重新拟合。脚本以原系数复算13点样本内R²并与既有指标核对。
- 未运行 CatMAP，未更改 TOF、DFT/NEB、反应网络或旧图。
- 补全方式：仅移除原图的凸包遮罩，在完整矩形上直接评估同一个二次多项式。矩形外推区只是经验模型数学延伸，不代表新催化剂结果。特别是若最大值落在扩展矩形边界且远离催化剂凸包，不得解释为物理活性峰。
- 形成能为历史参考下的 O、Ha、Hb 描述符；Ha 与 Hb 分开。历史单物种形成能与当前 CatMAP 能量参考尚未完成物理一致性验证。
- 新图标题明确说明为经验二次响应面、不是 CatMAP 连续网格。点为13个真实已验证催化剂；模型表面不是新的 CatMAP 计算。
- 网格 XYZ 每种分辨率均为 N×N 行且无 NaN。矩阵 CSV 左上角为 `x\y`，第一列为 x，首行其余列为 y，内部为预测 log10(TOF)。

## 坐标范围与最大值

| 图 | xmin (eV) | xmax (eV) | ymin (eV) | ymax (eV) | full-range maximum (x, y, log10TOF) | convex-hull maximum (x, y, log10TOF) |
|---|---:|---:|---:|---:|---|---|
| O*+Ha* | -3.5815130250 | 2.2445055250 | -4.0330441750 | 0.6444456750 | (-3.5815130250, -4.0330441750, 42.7173694816) | (1.9771240938, -2.0075455970, 2.0308891707) |
- O*+Ha* full-range maximum: boundary maximum; convex-hull maximum: boundary maximum. 两者不重合。全域最大值附近催化剂（标准化描述符距离最近）：Ti-Ni, Ti-Mo, Ti-Cr。
| O*+Hb* | -3.5815130250 | 2.2445055250 | -4.7225265250 | 0.7313290250 | (2.2445055250, -3.5702686858, 4.2497395310) | (1.9460959895, -3.5146867243, 2.4122719010) |
- O*+Hb* full-range maximum: boundary maximum; convex-hull maximum: boundary maximum. 两者不重合。全域最大值附近催化剂（标准化描述符距离最近）：Ti, Ti-Ni, Ti-Co。

## Origin 文件
- 精细绘图优先导入 `*_empirical_fullgrid_400x400_xyz.csv`：XYZ 长表最直接，三列为 x、y、predicted_log10_TOF。
- 若使用 Origin Matrix/Contour，可导入 `*_empirical_fullgrid_400x400_matrix.csv`；200×200版本用于较轻量绘图。
- 催化剂叠加点与峰值分别见 `*_catalyst_points.csv` 和 `*_peak_summary.csv`。峰值在连续解析二次面/凸包上的精确优化位置，不一定恰好落在离散网格节点。

## 输出文件
- `O_Ha_TOF_empirical_map_full.pdf`
- `O_Ha_TOF_empirical_map_full.png`
- `O_Ha_catalyst_points.csv`
- `O_Ha_empirical_fullgrid_200x200_matrix.csv`
- `O_Ha_empirical_fullgrid_200x200_xyz.csv`
- `O_Ha_empirical_fullgrid_400x400_matrix.csv`
- `O_Ha_empirical_fullgrid_400x400_xyz.csv`
- `O_Ha_peak_summary.csv`
- `O_Hb_TOF_empirical_map_full.pdf`
- `O_Hb_TOF_empirical_map_full.png`
- `O_Hb_catalyst_points.csv`
- `O_Hb_empirical_fullgrid_200x200_matrix.csv`
- `O_Hb_empirical_fullgrid_200x200_xyz.csv`
- `O_Hb_empirical_fullgrid_400x400_matrix.csv`
- `O_Hb_empirical_fullgrid_400x400_xyz.csv`
- `O_Hb_peak_summary.csv`

输入 matched_data.csv SHA256：`79638b616d58fd4fbf2c2b06b30752d4d4e1dbe58b19bd99ea38788d12a7d4a0`
输入 fit_metrics.json SHA256：`f899db02419088590d51f939372be361cd257395136e7ec2a7e9d4213ab824a2`
