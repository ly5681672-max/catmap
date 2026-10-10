# lyq_alkene_epoxidation1 — 独立 NEB 共吸附 CatMAP 模型

**独立目录，原版 `lyq_alkene_epoxidation` 原封不动。当前在 `feature/lyq-alkene-epoxidation1-neb`，不合并 main。**

## 文件及运行
- `己烯环氧化.xlsx`：独立拷贝的原始计算数据；读取 Sheet「吉布斯自由能汇总」。
- `energies.txt`：已生成的 15 个表面的完整共吸附态 **形成自由能** 和四气相物种参考值。无需重新读取 Excel 即可运行。
- `coadsorption_model.mkm`：完整共吸附反应网络的 CatMAP 输入。
- `run_model.py`：一键检查和调用 CatMAP；显式确认临时势垒才可求解。
- `build_coadsorption.py`：可选，从**本目录 Excel** 重建输入（需要 `artifact_tool`）。
- `energy_audit.csv`：原始 G 值、三态差值、暂定有效动力学势垒。
- `check_neb.py` 与 `neb_three_state_summary.csv`：最初与八步模型的辅助对照（前者跨目录对照，**不是新版运行依赖**）。
- `coadsorption_model.mkm.template`：历史初稿，仅供比较，不用于运行。

先激活安装了 CatMAP 的 Python 环境（如本地已有 `H:\\catmap\\catmap`，请自行确保模块在该环境可导入），在本目录运行：

```powershell
python run_model.py --validate-only
python run_model.py --allow-provisional
```

第二条会运行 CatMAP，并在本目录生成求解缓存/结果。如需重新由 Excel 提取（仅在已安装 artifact_tool 时）：

```powershell
python run_model.py --rebuild --validate-only
```

## 反应状态与元素守恒
- `R*`：共吸附 `(H2O2)(C6H12)*`，组成 C6H14O2。
- `P*`：共吸附 `(C6H12O)(H2O)*`，组成 C6H14O2。
- `RP*`：**暂定有效过渡态**；仅用于让 CatMAP 在非负势垒模型下工作，不是已验证的一阶鞍点。

```text
* + H2O2_g + C6H12_g <-> R*
R* <-> RP* -> P*
P* <-> C6H12O_g + H2O_g + *
```

相比旧版八步机制，这三步直接描述总体协同环氧化；不能直接将其动力学参数解释为旧版 OOH 或 H 转移基元反应的参数。

## 形成自由能参考
对于完整共吸附态 X∈{R,P,RP}：

```text
Gf(X*) = G_sheet(X*) - G_sheet(clean_slab) - 6*mu_C - 14*mu_H - 2*mu_O
Gf(gas) = G_sheet(gas) - sum(n_element * mu_element)
mu_C=-9.28, mu_H=-1.11, mu_O=-4.37 eV
```

严格说这是指定元素参考态上的形成自由能，不是传统热力学标准生成 Gibbs 能。元素基准只用于一致化，同原子计量反应的差值相消。选择冻结热力学项 `frozen_gas`/`frozen_adsorbate` 避免重复添加自由能校正。

读取行：裸板 4；R 12；路径标记 16；P 20（TiMo P = 第21行组合总能）；气相 G = 57–60；TiTa 的 G 列特殊偏移。对所有模型点按相同参考计算。

## 势垒及物理限制
`G_RP = max(G_R, G_marker, G_P)`（同一表面的绝对自由能）。
这是**用于数值求解的暂定势垒假设**。尤其 TiNb/TiW/TiTa 三态下降，不等于整条 NEB 图像经审查完全无峰；**只有 NEB 全图像 + 适当自由能修正验证后才可确认无能垒**。其余表面的单个 marker 也不等于真实最高鞍点。

气相后缀代表 CatMAP 的外部化学势库；液相 H2O2/C6H12/H2O/C6H12O 的 `pressure` 是**近似有效活度**，并非真实气相分压。模型使用单总位点、无显式溶剂作用及单位活度初期环氧产物近似。绝对 TOF 未经过实验/完整 NEB 验证，不应发表为定量预测。

**独立性**：运行/生成只使用本目录的 Excel、energies、配置和脚本；无需原有 lyq_alkene_epoxidation 文件。唯一外部软件依赖为 Python+CatMAP，选择重建时另需 artifact_tool。

## 状态
- 数据及源文件自包含，已通过静态能量守恒和形成能一致性检查。
- 已提供可直接调用 CatMAP 的运行入口。
- CatMAP 运行和收敛**需要在有 CatMAP 依赖的用户环境中执行，当前交付环境未安装 CatMAP，不能宣称已跑通求解器**。
