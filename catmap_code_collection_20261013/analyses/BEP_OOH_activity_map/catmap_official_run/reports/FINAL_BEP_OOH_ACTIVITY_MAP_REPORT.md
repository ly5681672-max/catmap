# CatMAP 0.3.1 四参数 OOH* 扩展BEP：完整计算结果

## 执行口径

- 本机CatMAP：0.3.1；实际调用 `ReactionModel.run(recalculate=True)`，模型配置为 `ThermodynamicScaler` + `SteadyStateSolver`，不修改CatMAP安装源码。
- 13个基线逐位点使用归档严格 `.mkm`/`energies.txt`的隔离快照；333.15 K、单有效中心、原六步反应、effective reservoir concentration（H2O2=0.185、C6H12=1.0、H2O=0.815、C6H12O=1e-20）、frozen气/吸附热力学、mpmath 75位、残差容差1e-70。
- OOH输入为历史单吸附态形成能，仅作外部BEP描述符；不等同六步网络桥连OOH，也没有加进网络。
- exact fit `R²=0.886665234`, RMSE=0.157751613 eV, MAE=0.121556863 eV, LOOCV RMSE=0.381431870 eV。系数/中心见 `data/exact_BEP_fit.json`；LOOCV每折只用12个位点并重算中心。

## CatMAP执行计数与验收

- 原始基线：13个位点；归档TOF最大绝对差 log10=2.2508e-11 dex。
- BEP全样本拟合重算：13个位点；LOOCV逐点重算：13个位点。
- 网格61×61：3721个坐标，CatMAP求解状态计数 `{"CONVERGED": 3721, "CONVERGED_ADAPTIVE": 0, "REUSED_CACHE": 0, "FAILED": 0, "FAILED_VALIDATION": 0}`；结果来源阶段计数 `{"exact_hash_cache_from_grid_5x5": 25, "validated_CatMAP_result_from_original_61x61_stage": 3696}`；成功=3721；最大残差=8.422906e-71；容差=1e-70。其中来自61×61正式阶段的3696个节点均有该阶段逐节点记录；其余与5×5节点坐标及输入哈希完全一致的25个节点复用先前已验证CatMAP结果。当前报告重导出时3721行均由精确输入哈希缓存读取，见`retrieval_status`；这不改变`status`所记录的真实求解结果。
- 5×5和15×15阶段数据分别保存为 `data/grid_5x5.csv`、`data/grid_15x15.csv`。所有失败坐标不填零/不插值，日志见 `runs/`。
- 61×61正式网格求解阶段耗时=337.11 s（原始计算调用计时，不含后续缓存复用的图件重生成）；本次报告/图件重生成调用耗时=21.33 s。此前输入准备、基线、BEP和5×5/15×15阶段为独立调用，未合并计入该值。全域最大值分类=`Maximum within sampled range (boundary)`；TS1高于两个端点的子域最大值分类=`Maximum within independent-TS-valid sampled range (boundary)`。
- TS1高于初/末态的网格点=2779/3721；端点控制节点不删去，图中以斜线标识。
- 每个节点均使用运行前绑定的 `catmap_case.mkm` 指向本节点 `energies.txt`；解析后逐项比对Rpair/TS1/EpoxideHHO/HHO/TS2/OHHa/H2O共7态。本工作流修复前生成的旧哈希运行/缓存已保留作诊断，但不进入当前CSV、峰值或图件；`RUN_SCHEMA`变化使其不能作为当前结果缓存。

## 网格最大值坐标

- 全部网格最大值：ΔG₂=-0.179686 eV，E(OOH*)=-6.232301 eV，log10TOF=1.476181907191137478389875513320083050134310222809013260820659319992020060029044850064597578873048788；状态=`Maximum within sampled range (boundary)`，TS1独立过渡态有效=`False`。
- TS1高于初态和末态的子域最大值：ΔG₂=-0.179686 eV，E(OOH*)=-5.880416149999999 eV，log10TOF=1.476181907184934344689900662219428826086354846811761875186798467377142175407431293446049507289128727；状态=`Maximum within independent-TS-valid sampled range (boundary)`。
- 两者均为边界最大值；本网格没有证据支持把最大值称为内部火山峰。

## TS1物理有效性与CatMAP端点处理

保留所有BE预测的有符号 TS1−Rpair 和预测TS1原值。独立过渡态有效条件是 `Y1 >= max(0, ΔG2)`。本机 CatMAP 0.3.1官方 `catmap/data/templates.py` 的 `rate_constants` 实现以 `max(G_IS, G_FS, G_TS)` 作为有效势垒顶，因此当预测TS低于任一端点，实际速率由较高端点控制；这不是输入能量被悄悄改写。本结果分别保留 `independent_TS_above_both_endpoints`、`CatMAP_endpoint_clamp_applied` 和 `CatMAP_effective_TS_control`。

## 条件网格限制与科学解释

网格保持Ti-Co的Rpair、HHO、OHHa、H2O、TS2和气相储库固定，只随ΔG2改变EpoxideHHO，并用历史单吸附态OOH值经四参数扩展BEP预测TS1。由此得到的是Ti-Co条件背景下的CatMAP网格，不是13种真实催化剂分别重算的排名。OOH历史参考与严格六步形成能的slab/热修正同源性仍未完成物理闭合；拟合R²不能替代可迁移验证，LOOCV误差已并列给出。BEP点位重建造成的TOF变化及网格峰，均不能单独证明双核协同增强或普适火山规律。

## 关键输出

- `data/baseline_validation.csv`、`data/13site_original_TOF.csv`、`data/13site_BEP_fitted_TOF.csv`、`data/13site_BEP_LOOCV_TOF.csv`
- `data/BEP_TS1_predictions_and_LOOCV.csv`、`data/grid_61x61.csv`、`data/grid_61x61_coverages.csv`、`data/grid_61x61_net_rates.csv`
- `figures/`下3类PNG/PDF/SVG图；输入SHA256见 `data/SHA256SUMS.csv`。

## Ti和Ti-Co实际TOF

|位点|原始CatMAP TOF|全样本BEP TS1 CatMAP TOF|LOOCV-BEP TS1 CatMAP TOF|
|---|---:|---:|---:|
|Ti|7.691724096364e+00|9.253338370938e-01|5.001015383036e-01|
|Ti–Co|3.708424176046e+00|3.677386351576e+00|3.804940306853e-13|
