from pathlib import Path
import json, hashlib
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
D,O,F,R=ROOT/"data",ROOT/"origin_exports",ROOT/"figures",ROOT/"reports"
S=json.loads((D/"bep_analysis_summary.json").read_text(encoding="utf-8"))
MET=pd.read_csv(O/"official_BEP_comparison.csv")
EXT=pd.read_csv(O/"BEP_extended_model_ranking.csv")
MAP=pd.read_csv(O/"BEP_TOF_mapping_summary.csv")
RES=pd.read_csv(O/"BEP_residual_analysis.csv")
cat_meta=json.loads((D/"official_BEP_scaler_metadata.json").read_text(encoding="utf-8"))

def row(path): return MET.loc[MET.path.eq(path)].iloc[0]
f1=row("NEB1_path2"); f2=row("NEB2_path4")
def fmt(v,n=4): return f"{float(v):.{n}f}"
neg1=RES.loc[(RES.path.eq("NEB1_path2"))&(RES.signed_TS_minus_IS_eV<0),["catalyst","signed_TS_minus_IS_eV"]]
neg2=RES.loc[(RES.path.eq("NEB2_path4"))&(RES.signed_TS_minus_IS_eV<0),["catalyst","signed_TS_minus_IS_eV"]]
ext1=EXT.loc[EXT.path.eq("NEB1_path2")].sort_values("RMSE_eV").iloc[0]
ext2=EXT.loc[EXT.path.eq("NEB2_path4")].sort_values("RMSE_eV").iloc[0]
within=MAP.loc[MAP.existing_grid_bounds_supported]

(R/"BEP_EXISTING_WORK_AUDIT.md").write_text(f"""# 既有BEP工作与数据入口审计

## 本轮复用的数据

- 严格CatMAP形成能表快照：`data/source_strict_energies.txt`，按`surface_name`与13个`CONVERGED`位点匹配；Ti–Zr排除。
- 13位点TOF及历史O/Ha/Hb：`data/source_matched_O_H_TOF.csv`。
- 既有BEP逐点值、回归摘要和拟合链预测：`data/source_prior_BEP_points.csv`、`data/source_prior_BEP_summary.csv`、`data/source_prior_chain_predictions.csv`。
- 含碳/单吸附态候选及参考说明：`data/source_carbon_O_H_TS1.csv`、`data/source_descriptor_inventory.csv`。
- 已有CatMAP 50×50直接求解网格和400×400 VectorMap插值：`data/source_catmap_grid_50x50.csv`、`data/source_catmap_VectorMap_400x400.csv`，连同`data/source_grid_summary.json`。
- 官方Scaler隔离核验产物：`data/official_BEP_scaler_coefficients.csv`与`data/official_BEP_scaler_metadata.json`。

## 数据流核验

13个位点由`surface_name`键连接，不用行号拼接。BEP输入由严格表中反应物/产物/TS形成能重新计算，并与既有逐点表核对。NEB有符号TS−IS差保留负值。历史单吸附态数据只作为探索性扩展特征，参考闭合未验证，不能解释成已统一到严格模型slab。

## 范围与禁止操作

本包不调用`ReactionModel.run()`，不重跑稳态、不扫条件、不更改DFT、NEB、TOF、反应机理或历史网格。CatMAP连续面仅复用已存VectorMap插值；不是新求解的CatMAP网格。本轮代码与输出先在C盘隔离目录生成，完成检查后再作为新目录复制入H盘归档；原始数据不覆盖。
""",encoding="utf-8")

(R/"BEP_SCALER_VERIFICATION.md").write_text(f"""# CatMAP 0.3.1 官方BEP Scaler隔离核验

- 本机版本：{cat_meta.get('catmap_version','0.3.1')}。
- 使用官方`GeneralizedLinearScaler`的BEP转换矩阵路径；调用`get_coefficient_matrix()`以让官方Scaler按正常调用链初始化约束和矩阵状态。
- 只参数化TS标度关系；元数据确认`ReactionModel.run()`未调用，未生成稳态或TOF。
- TS关系采用初态基准：`E_TS − E_IS = α(E_FS − E_IS) + β`。本项目拟合使用相同有符号定义。

| 路径 | CatMAP α | NumPy OLS α | CatMAP β (eV) | NumPy OLS β (eV) | 样本内 R² | RMSE (eV) | LOOCV RMSE (eV) |
|---|---:|---:|---:|---:|---:|---:|---:|
| NEB1 / path 2 | {fmt(f1.CatMAP_alpha,8)} | {fmt(f1.numpy_alpha,8)} | {fmt(f1.CatMAP_beta,8)} | {fmt(f1.numpy_beta,8)} | {fmt(f1.R2)} | {fmt(f1.RMSE_eV)} | {fmt(f1.LOOCV_RMSE_eV)} |
| NEB2 / path 4 | {fmt(f2.CatMAP_alpha,8)} | {fmt(f2.numpy_alpha,8)} | {fmt(f2.CatMAP_beta,8)} | {fmt(f2.numpy_beta,8)} | {fmt(f2.R2)} | {fmt(f2.RMSE_eV)} | {fmt(f2.LOOCV_RMSE_eV)} |

CatMAP和独立OLS系数在浮点精度内一致。隔离模型配置和调用标志见JSON；模型设置与真实严格运行条件不相混，且没有执行动力学。
""",encoding="utf-8")

(R/"BEP_MECHANISM_REPORT.md").write_text(f"""# BEP关系、共吸附修正与TS重建：结果报告

## 结论摘要

| 路径 | BEP关系（有符号） | R² | RMSE (eV) | LOOCV RMSE (eV) | LOOCV R² |
|---|---|---:|---:|---:|---:|
| NEB1 / Rpair→TS1→EpoxideHHO | `ΔE‡signed = {fmt(f1.numpy_alpha,8)} ΔG_rxn + {fmt(f1.numpy_beta,8)}` | {fmt(f1.R2)} | {fmt(f1.RMSE_eV)} | {fmt(f1.LOOCV_RMSE_eV)} | {fmt(f1.LOOCV_R2)} |
| NEB2 / HHO→TS2→OHHa | `ΔE‡signed = {fmt(f2.numpy_alpha,8)} ΔG_rxn + {fmt(f2.numpy_beta,8)}` | {fmt(f2.R2)} | {fmt(f2.RMSE_eV)} | {fmt(f2.LOOCV_RMSE_eV)} | {fmt(f2.LOOCV_R2)} |

`ΔG_rxn = G_FS − G_IS`，`ΔE‡signed = G_TS − G_IS`。这两个关系是13个既有位点上的样本内拟合；LOOCV表现较弱，尚不能作为可靠外推预测器。NEB1负的原始TS−IS值出现在：{", ".join(f"{r.catalyst} ({r.signed_TS_minus_IS_eV:.6f} eV)" for r in neg1.itertuples())}。NEB2负值出现在：{", ".join(f"{r.catalyst} ({r.signed_TS_minus_IS_eV:.6f} eV)" for r in neg2.itertuples())}。这些是输入数据的有符号端点差，不被替换成零，也不被重新命名为正活化能。

## 单额外共变量的诊断

在基准BEP反应能之外逐次加入一个候选能量特征。当前样本内RMSE最低的候选：NEB1为`{ext1.added_descriptor}`（N={int(ext1.N)}，R²={fmt(ext1.R2)}，RMSE={fmt(ext1.RMSE_eV)} eV，LOOCV RMSE={fmt(ext1.LOOCV_RMSE_eV)} eV）；NEB2为`{ext2.added_descriptor}`（N={int(ext2.N)}，R²={fmt(ext2.R2)}，RMSE={fmt(ext2.RMSE_eV)} eV，LOOCV RMSE={fmt(ext2.LOOCV_RMSE_eV)} eV）。该候选为历史单吸附态OH*能量，当前未与严格模型slab统一；故仅是探索性相关，不据此宣布物理修正已成立。对所有特征的定义与参考状态限制见`origin_exports/BEP_extended_model_ranking.csv`。NEB1首步的`ΔG1=G(Rpair*)−G(H2O2_g)−G(C6H12_g)`是整体共吸附形成热力学差，不是纯粹的共吸附相互作用能。

## Ti–Fe敏感性与共吸附解释

去掉Ti–Fe后，NEB1基准R²={fmt(S['NEB1_path2_TiFe_excluded_sensitivity']['R2'])}、RMSE={fmt(S['NEB1_path2_TiFe_excluded_sensitivity']['RMSE_eV'])} eV；NEB2基准R²={fmt(S['NEB2_path4_TiFe_excluded_sensitivity']['R2'])}、RMSE={fmt(S['NEB2_path4_TiFe_excluded_sensitivity']['RMSE_eV'])} eV。该删除敏感性检查仅作诊断，正式主结果仍包括全部13个位点。严格Rpair态相对于气相反应物的ΔG1已导出；其解释包含各自吸附与共吸附整体效应，不能独立拆成“相互作用能”。

## TS重建与既有CatMAP表面映射

`BEP_TS1_TS2_reconstruction.csv`分别保留使用全样本拟合和LOOCV barrier预测、再加回真实对应IS能量的TS重建结果。前者是in-sample诊断，后者是更保守的留一预测。另复用既有O/H描述符链的LOOCV预测（NEB1 Ha、NEB2 Hb），按其TS1/TS2坐标查询既有400×400 VectorMap插值表：{len(within)}/{len(MAP)}个位点坐标落在保存网格范围内。导出的TOF仅为最近插值网格点值，不是新CatMAP求解，也不是重建TS的真实TOF；网格外不赋值。底图的极值或网格质量仍以原有50×50直接求解审计为准。

## 科学判断与限制

1. 两条传统BEP线性关系在样本内解释度偏低，LOOCV不支持稳健跨位点预测。不能把全局线性BEP当作已验证普适定律。
2. NEB1样本内拟合误差明显较大；NEB2绝对误差较小但LOOCV仍有限。负的TS−IS端点差保留，可能反映局部端点能量/构型排序，不应被端点差直接解释成通常意义的正活化势垒。
3. OH*额外项在两条反应的样本内拟合均有改善。NEB1的LOOCV RMSE由{fmt(f1.LOOCV_RMSE_eV)}增至{fmt(ext1.LOOCV_RMSE_eV)} eV，未改善；NEB2则由{fmt(f2.LOOCV_RMSE_eV)}降至{fmt(ext2.LOOCV_RMSE_eV)} eV，呈现有限的留一改善。但OH*历史slab参考闭合不完整，样本仅13个位点，因此尚不能升级为正式共吸附修正或外部有效的BEP模型。
4. 旧O/H单物种数据与严格六步态能量存在位点相关参考偏移；尤其Ti–Fe历史clean-slab差异约1.064927 eV（历史审计记录），当前缺少逐原子结构、cell和同一热修正口径的证据来解释该差值。因此旧吸附能仅作探索性变量，不能与严格Rpair/TS绝对能量直接混成统一参考。
5. 首步共吸附状态的真实结构与独立吸附态简单叠加并不等价。模型中的Rpair自由能可以作为整体状态量，但当前资料不足以单独分离吸附能和吸附物间相互作用能。
6. 官方CatMAP B.E.P. scaler系数已在隔离模式下复现；它验证的是实现/拟合系数，不验证拟合质量或物理普适性。CatMAP文档所述BEP是TS能相对初态对反应能的标度关系，具体解释仍依赖一致的参考能和物种定义。

## 复现与产物

- `scripts/verify_catmap_bep.py`：仅调用官方Scaler系数矩阵，无`ReactionModel.run()`。
- `scripts/analyze_bep_reconstruction.py`：重建13位点BEP、扩展项、留一评估、TS重建和既有VectorMap近邻读取。
- `scripts/render_bep_figures.py`：科研图及Origin表输出。
- `data/official_BEP_scaler_metadata.json`：注明实际调用状态。
- `reports/SHA256SUMS.csv`：本归档包文件哈希。

未做：新的DFT/NEB、CatMAP稳态/TOF、反应条件搜索、CatMAP网格计算。
""",encoding="utf-8")

files=[]
for p in sorted(ROOT.rglob("*")):
    if p.is_file() and p.name!="SHA256SUMS.csv":
        h=hashlib.sha256(p.read_bytes()).hexdigest()
        files.append({"relative_path":str(p.relative_to(ROOT)),"sha256":h,"bytes":p.stat().st_size})
pd.DataFrame(files).to_csv(R/"SHA256SUMS.csv",index=False)
print("Reports written; package files hashed:",len(files))
