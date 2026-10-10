# NEB 共吸附独立对照版说明

与原目录的子目录、文件名称和相对位置保持相同。**此副本**对应整体共吸附反应；旧项目仍为八步基元机理。二者的 `energies.txt` 不能混用。

## 运行

1. 激活可导入 `catmap` 的 Python 环境。
2. `cd lyq_alkene_epoxidation1/2-creating_microkinetic_model`。
3. `python -c "from catmap import ReactionModel; ReactionModel(setup_file='alkene_epoxidation.mkm').run()"`。
4. 在 `test.ipynb` 中复现和检查。

## 生成数据

在 `1-generating_input_file` 执行 `python generate_input1.py`；其输入直接来自本目录上一层的 Excel，写入同目录的 `energy_audit.csv` 及模型目录的 `energies.txt`。额外依赖 `artifact_tool`，如未安装可直接使用已提交的能量表。

## 反应与物种

- `R*`: (H2O2)(C6H12) 完整共吸附态，C6H14O2。
- `P*`: (C6H12O)(H2O) 完整共吸附态，C6H14O2。
- `RP*`: 按 `max(G_IS,G_marker,G_FS)` 建立的**数值有效 TS**；未验证为 CI-NEB 最高点。
- `* + H2O2_g + C6H12_g <-> R*`
- `R* <-> RP* -> P*`
- `P* <-> C6H12O_g + H2O_g + *`

气相后缀的 pressure 参数在此仅为液相近似有效活度；TOF 尚未通过 CatMAP 实跑及完整 NEB 校准。不同自由能修正基准不可混用。TiMo 的末态取 Excel 第21行组合总能，TiTa 的 G 列单独偏移。

## 与原版直接比对

- 原：`lyq_alkene_epoxidation/2-creating_microkinetic_model/alkene_epoxidation.mkm`
- 新：`lyq_alkene_epoxidation1/2-creating_microkinetic_model/alkene_epoxidation.mkm`

项目根目录及子目录命名仅差一个 `1`，其他主要对应文件路径保持一致。

## 2026-10-10 运行日志排查与绘图修复

- 已检查运行日志：`mapper_iteration_1: status - 0 points do not have valid solution`，表明本次 10×10 网格最终获得数值解，初始失败点经重试/二分恢复；这不等于结果已物理验证。
- `Damping step size exceeded max_damping` 与 `stationary point or singular jacobian` 是数值刚性/收敛尝试告警。本次没有盲目调整 DFT 能量或把它们隐藏。
- `header_evaluation: fail - could not save ...` 来自 CatMAP 生成文本日志时无法执行某些 NumPy/mpmath 表示。已在 Notebook 第一单元使用 CatMAP 自身的 `_pickle_attrs` 将这些属性写入 .pkl，避免同类文本重建失败；这不改变速率常数的计算。
- 原版所具有的速率、产物生成速率、覆盖度和 scaling PDF 图，已在新版 `2-creating_microkinetic_model/test.ipynb` 中逐图恢复。Notebook 增加 100 点网格完整性检查、覆盖度值域检查、结果表及 `ts_scaling_audit.csv`。
- **核心物理差异**：`RP_s` 使用 GeneralizedLinearScaler 的线性拟合，模型并不逐催化剂保留 `energies.txt` 输入的 `RP`。在用户实际日志中 `G_RP^{fit}≈0.87204633·G_R−4.01396471`。据此，TiW/TiNb/TiTa 的输入 0 eV 有效势垒，在标度映射中分别变为约 0.183/0.242/0.259 eV；该拟合误差对 TOF 可造成显著影响。`ts_scaling_audit.csv` 显示所有 15 个催化剂的差别，**切勿把其 TOF 当作精确 NEB 输入后的动力学预测**。
- 修复了 `1-generating_input_file/generate_input1.py` 的 Excel 路径错误，使其只读取 `../己烯环氧化.xlsx`。
- 本仓库的更改已通过文件核对、工作表数值一致性与 Notebook 静态语法检查，**尚未在用户本地 CatMAP 环境完成绘图和结果收敛的动态复核**。

### 本地重新运行

在 `lyq_alkene_epoxidation1/2-creating_microkinetic_model` 打开 `test.ipynb`，按单元 1→7 顺序执行。不要使用旧版 `lyq_alkene_epoxidation` 的工作目录。预期新增 `rate.pdf`、`all_rates.pdf`、`production_rate.pdf`、`pretty_production_rate.pdf`、`coverage.pdf`、`R_P_coverage.pdf`、`scaling.pdf`、`production_rate_table.txt`、`ts_scaling_audit.csv`。第 8 单元是默认注释的敏感性计算选项，以免在基础验算时增加高刚性求解工作。

## 2026-10-10 最新实跑错误（已修改 Notebook，待用户本地复跑）

GitHub 最新已执行的 `test.ipynb`（In[12]–In[19]）包含以下可复现问题：

1. **In[13]：TypeError**：`model.resolution` 是列表 `[10,10]`，误用了 `int(model.resolution)`。现改为对每维取乘积，预期 100 个描述符点。
2. **In[15]：PermissionError**：`rate.pdf` 已存在且可能被 PDF 预览器锁定，Matplotlib 无法以 `wb` 覆盖。新版绘图首先尝试原文件名，受锁时回退到 `plots/rerun_YYYYmmdd_HHMMSS/<原文件名>`；若目录不可写仍需排查权限。
3. **In[15]：RuntimeWarning invalid value encountered in log10**：净基元速率有符号，不能直接按对数绘制。新版 `rate.pdf`、`all_rates.pdf` 改为线性带符号速率；环氧产物净生成速率只有在全部非负时才采用对数显示。
4. **CatMAP 文本日志序列化**：`coverage_map` 等大型映射数据无法被 CatMAP 的 `exec(repr(...))` 直接解析，新增 CatMAP 原生 `_pickle_attrs` 存档设置，以使它们从 .pkl 重建。旧日志仍保留以前的警告，须在重跑后复核。
5. **实际已得到结果**：最新 `production_rate_table.txt` 已有 100 个网格点和 4 种气相物种的生成速率；`ts_scaling_audit.csv` 共 15 组催化剂，11 组 RP 拟合偏差 >0.05 eV。零有效输入势垒的 TiNb/TiW/TiTa 在当前拟合中对应 0.2416/0.1829/0.2586 eV。此问题属于动力学映射假设，**目前尚未通过编程修正为真实逐催化剂 NEB 势垒**。

### 修复后的本地复跑顺序

```powershell
cd H:\\catmap\\catmap
git status --short
git pull --ff-only origin main
```

确认 PyCharm/Jupyter 刷新新版 `lyq_alkene_epoxidation1/2-creating_microkinetic_model/test.ipynb`，选择重新启动 Kernel、顺序运行各代码单元。首先关注网格完整性结果 `production_rate=100, rate=100, coverage=100` 和物理覆盖度警告数量；其次检查所有 PDF 保存提示（遇到锁定时自动转 `plots`），最后检查 TS 拟合误差文件。若缓存很大，强制重新计算可能花费一些时间。
