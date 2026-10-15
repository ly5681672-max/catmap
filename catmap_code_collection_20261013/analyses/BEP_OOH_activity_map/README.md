# 基于四参数OOH*扩展BEP的二维活性图（Ti–Co条件参照）

## 结果与口径
本目录提供新的**ΔG₂–E(OOH*)二维活性预测**，不再沿用原TS1–TS2双描述符火山图。BEP关系为：

```
u = ΔG₂ + 2.0922766923 eV
z = E(OOH*) + 4.3098743846 eV
Y1_pred = G(TS1*) − G(Rpair*) =
  0.4121230973 + 0.3703773006*u + 0.2002174961*z − 0.3097773970*z²
```

13位点NEB1有符号TS−IS能量差样本内BEP拟合R²=0.886665；原始OOH*为**历史单吸附态形成能**，不能直接称为与严格6步模型完全同源的桥连OOH描述符。

**明确的条件模型：** 以Ti–Co为参照，固定Rpair*、HHO*、OHHa*、H₂O*、TS2的能量及全部气相储库；网格中改变\(\Delta G_2\)时只改变\(G(\mathrm{EpoxideHHO}^*)=G(\mathrm{Rpair}^*)+\Delta G_2\)，同时通过扩展BEP重建TS1。对于预估TS1低于初始态的网格，使用CatMAP式端点夹断以维持非负势垒，灰色显示、不得当作真实NEB过渡态。其他12个位点标注**原始描述符坐标**，不是本条件面在不同催化剂的预测TOF。

**执行口径：** 图背景采用对CatMAP六步单有效活性中心稳态方程的独立解析重求解（树定理/循环守恒，采用原归档的能量、储库浓度和Eyring速率常数），**不是**新的官方CatMAP 0.3.1程序网格运行。不能用这张图声称双核已超过单Ti或出现物理上可靠的火山峰。独立解析模型先以原Ti–Co的真实TS1能量重现归档速率，误差<1e-4 log10(TOF)。

## 基准与网格数据
- Ti–Co档案log10TOF：0.569189404；独立基线复现：0.569188721。
- 将Ti–Co的TS1替换成新BEP重建能量后，在原ΔG₂/OOH坐标下得到log10TOF：0.565538591，这是条件预测而非新官方CatMAP计算。
- 61×61=3721个独立解析节点；端点夹断比例：0.253。
- 完整网格数值最高点：ΔG₂=-0.179686, EOOH=-6.232301, log10TOF=1.476181；未夹断区最高点：ΔG₂=-0.179686, EOOH=-5.880416, log10TOF=1.476181。最大值位于网格边界，**不能视为发现最优催化剂或真实火山峰**。

## 文件
- `figures/01_BEP_OOH_dG2_TiCo_conditional_activity.svg`：SVG科研图（可用于PPT）。
- `data/grid_BEP_OOH_dG2_TiCo_conditional.csv`：全部网格数值及是否夹断。
- `scripts/plot_from_grid.py`：使用归档CSV重绘PNG、PDF、SVG。
- `scripts/check_native_catmap_nodes.py`：准备并执行少量原生CatMAP节点复核（需用户本机CatMAP 0.3.1；**本次尚未运行**）。

运行重绘：
```powershell
cd H:\catmap\catmap\catmap_code_collection_20261013
python .\analyses\BEP_OOH_activity_map\scripts\plot_from_grid.py
```

这张图目前可用于说明“新BEP预测TS1与ΔG₂变化影响Ti–Co条件反应网络TOF”，不能替代下一阶段：统一历史OOH能量参考、进行官方CatMAP节点验证、分别对13个位点用BEP新TS1重求TOF、验证完整标度网络及稳健性。
