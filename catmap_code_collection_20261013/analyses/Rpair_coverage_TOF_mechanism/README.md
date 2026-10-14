# Rpair* 共吸附态—TOF关联及动力学敏感性（13位点）
 
## 目的
检验实验关注的反应前驱共吸附态 `Rpair*` 与环氧产物TOF的模型内关系，**区分稳态覆盖度关联和真实速率敏感性**。当前六步网络的Rpair*定义为`H2O2 + C6H12`共吸附态，并不能未经核实直接等同于实验检测的桥连Ti–OOH中间体。

## 数据来源与方法
- 严格CatMAP 0.3.1基准；333.15 K；13个`CONVERGED`催化位点，排除未收敛Ti–Zr。固定储库：H₂O₂=0.185、烯烃=1、H₂O=0.815、环氧己烷=1e-20；这些量在原`.mkm`内以浓度输入，**不应擅自改称bar或mol/L**。
- `../microkinetic_mechanism_13sites/data/coverage_13sites.csv`：Rpair*稳态覆盖度。
- `../microkinetic_mechanism_13sites/data/rate_control_matrix.csv`：CatMAP官方产品TOF的DRC。
- `../microkinetic_mechanism_13sites/data/reaction_energy_13sites.csv`：ΔG₁。
- `../dG3_coverage_DRC_mechanism/data/combined_dG3_coverage_DRC_TOF.csv`：TOF和ΔG₃。
- 本目录`data/rpair_13sites_joined.csv`保存已配对数值（site ID一对一）。

## 已完成结果（来自已归档数据的后处理；不是新CatMAP求解）
- `figures/01_Rpair_coverage_vs_TOF.svg`：13位点log₁₀θ(Rpair*)—log₁₀TOF拟合，样本内R²=0.9365，Spearman ρ=0.9505，最小二乘斜率=0.8091，截距=5.0400。
- `figures/02_Rpair_DRC_heatmap.svg`：Rpair* / TS1* / EpoxideHHO* / HHO* / OHHa*的官方DRC，显示Rpair*数值近零。
- ΔG₁单变量拟合TOF：R²=0.2044；ΔG₃单变量拟合TOF：R²=0.9315（使用CatMAP有效ΔG₃，固定条件下与标准值仅差共同平移）。
- 覆盖度由**同一稳态模型**求得，故覆盖度—TOF高R²属于模型内相关，不能当作Rpair能量控制TOF的独立证据。

## 需要实际计算的下一阶段（仅提供脚本，尚未运行）
`scripts/run_rpair_sensitivity.py`读取`catmap_official_run/results/official_strict_tolerance/<site>`内的原版`.mkm`与`energies.txt`，在单独的输出目录复制输入，做少量Rpair*能量±0.025/±0.05 eV及H₂O₂、烯烃浓度的单因素扰动，每个条件调用**原生CatMAP 0.3.1**重新求稳态、TOF与覆盖度。脚本不改源模型。

在具有CatMAP 0.3.1环境的Windows PowerShell下：
```powershell
cd H:\catmap\catmap\catmap_code_collection_20261013
G:\anaconda\python.exe .\analyses\Rpair_coverage_TOF_mechanism\scripts\run_rpair_sensitivity.py --surfaces ti tico tife --python G:\anaconda\python.exe
```
输出`results/scan_rows.csv`及`results/energy_sensitivity.csv`（仅成功计算的数值），失败保留`FAILED`状态和日志。若基准计算与归档差异超阈值，自动中止；不把失败点当作零。

**解释提醒**：官方DRC很小，必须接受可能得出Rpair*并非当前条件下主要控制态的结果；单Ti/异核反应网络能量敏感性可能不同。能量扰动时保持TS1绝对能量固定，反应能和正反势垒会一同变化，符合状态能扰动导数概念，但不能把其等同于真实材料独立调节Rpair*的能力。
