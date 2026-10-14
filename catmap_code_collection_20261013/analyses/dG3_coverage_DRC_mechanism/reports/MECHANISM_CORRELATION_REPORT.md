# ΔG3–覆盖度–DRC–TOF及TS1–TS2峰值机制关联

## 数据与口径

使用CatMAP 0.3.1严格基准、333.15 K、13个已收敛位点；Ti–Zr未纳入。13位点数据按surface_name一对一合并。标准态ΔG3按形成能直接计算：

`ΔG3_standard = Gform(HHO*) + Gform(C6H12O_g) - Gform(EpoxideHHO*)`。

有效ΔG3取既有CatMAP固定储库严格基准的step-3结果；它包含气相化学势修正，与标准ΔG3的定义不同。两者在当前固定产物储库下相差约-1.322081875 eV且为常数，因此相关系数数值相同，不是两种独立证据。

## 13位点关联结果

- 标准ΔG3与log10TOF：Spearman ρ=-1.000，Pearson r=-0.965（N=13）。这复现了此前严格基准的单调相关，但本身不证明ΔG3是唯一控速因素。
- 标准ΔG3与EpoxideHHO覆盖度：Spearman ρ=0.893；与HHO覆盖度ρ=-0.668；与高精度空位覆盖度ρ=-0.502。
- EpoxideHHO覆盖度与log10TOF：Spearman ρ=-0.893。多个位点覆盖度极接近1，需结合完整数值和DRC，不把近饱和等同于唯一限速。
- 13位点中的多数双金属位点（不含Ti–Fe的11个）有11/11同时满足DRC(HHO)>0.9和DRC(EpoxideHHO)<-0.9；这支持局部活性对两种能量的相反敏感性。DRC仍是局部导数，不是唯一决速步骤标签。
- Ti–Fe不同：θ(OHHa*)≈0.910且DRC(OHHa*)≈-0.906，EpoxideHHO覆盖仅≈0.087。这与OHHa积累/后段敏感的特殊稳态响应一致；同时其clean-slab来源差异仍未物理闭合，故不把它提升为已验证的本征机制因果结论。
- Ti单核位点空位约0.998818，DRC(TS1)≈0.999，与多数双金属位点不同。

因此较稳妥的解释是：在多数异核位点上，较不利的环氧产物脱附端点热力学与强EpoxideHHO稳定/覆盖、HHO/EpoxideHHO相反的DRC共同形成活性趋势；Ti–Fe及Ti表明还有OHHa或TS1敏感的不同响应。第3步采用CatMAP端点势垒近似，ΔG3–TOF关系部分由该模型定义耦合，不能作为独立验证的普适描述符规律。

## TS1–TS2既有CatMAP网格峰及两侧

既有官方50×50 CatMAP缓存给出离散峰TS1=-33.314735 eV、TS2=-7.673499 eV、log10TOF=4.821900。本轮仅复用峰及其四个已缓存轴向邻点，并调用官方Scaler对这些坐标作能量读取；官方能量API在13个真实位点上与strict baseline最大差：ΔG 1.2e-07 eV、正向势垒 1.2e-07 eV。

- TS1轴：log10TOF范围 4.572956 至 4.821900；ΔG2/ΔG3/ΔG4与正向显式TS有效势垒见`data/TS1_TS2_peak_mechanism_data.csv`。
- TS2轴：log10TOF范围 4.335302 至 4.821900；ΔG2/ΔG3/ΔG4与正向显式TS有效势垒见`data/TS1_TS2_peak_mechanism_data.csv`。

峰处θ(EpoxideHHO)=0.123、θ(H2O)=0.477、θ(Rpair)=0.395；邻点显示TOF下降伴随这三种状态重新分配，而非单一覆盖物种替换可以充分概括。TS1两侧ΔG2/forward barrier和TS2两侧ΔG3/ΔG4变化均列在`data/TS1_TS2_peak_mechanism_data.csv`。

这些虚拟坐标的物种能量来自CatMAP标度器，不是新增DFT；也没有缓存该峰附近的DRC，因此峰两侧的机制解释是能量—覆盖度—TOF的探索性关联，不能声称已由局部DRC直接验证。

## 全体/12位点敏感性与限制

`data/TIFE_13site_vs_12site_comparison.csv`保留全体13位点，并标记包含Ti–Fe与否；12位点结果只作敏感性对照。排除Ti–Fe后标准ΔG3–log10TOF Spearman ρ=-1.000，与13位点趋势一致；这不构成排除Ti–Fe的理由。

图件包括13位点联合关系、覆盖度/DRC、峰值轴向覆盖度与TOF、峰值有效反应能，以及Ti–Fe代表位点对照；对应CSV均在`origin_exports/`。未重新绘制此前已有的单变量活动图，也未计算虚拟网格DRC。

## 科研表达分级

**可直接用于PPT（注明模型条件）：** 13位点严格基准的标准/有效ΔG3区分、覆盖度与TOF排序、官方DRC局部响应，以及CatMAP缓存峰/两侧的覆盖重分配；Ti–Fe具有不同的数值稳态响应。

**需补原始结构/计算证据后再作更强物理结论：** Ti–Fe历史与严格slab为何相差约1.065 eV；NEB2若干Gibbs修正的温度记录来源；峰区的DRC或更密的局部数据。上述缺口不影响本轮复用既有基线结果，但限制因果归因。
