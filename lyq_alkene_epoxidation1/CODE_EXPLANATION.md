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

## 2026-10-10 第二次实跑：绘图异常与 C/O/H 火山图

本次直接读取 GitHub 最新已执行 Notebook 输出：
- CatMAP 网格计算：100/100 点；coverage：0 个物理范围异常。
- Notebook rate/coverage 绘图：`TypeError: ufunc 'isfinite' not supported`（CatMAP 高精度 mpmath 结果引入 NumPy object 数组）。
- 原 `ScalingAnalysis`：`AttributeError: 'NoneType' object has no attribute 'keys'`（当前绘图模块接收到 `parameter_dict=None`）。
- `ts_scaling_audit.csv` 仍有 11/15 个催化剂 RP 拟合偏差 >0.05 eV。未修改 DFT 原始能量，未伪称拟合后的 TOF 为真实 NEB 值。

修复：
1. 在 `2-creating_microkinetic_model/test.ipynb` 中采用 NumPy `dtype=float64`，将 CatMAP `mpf` 标量转为可用于 Matplotlib 的数值。原 PDF 图名保持不变，输出包含 `volcano_R_P.pdf`（**仅 R/P 描述符**）。
2. 以 `RP_input` 对 `RP_scaled` 进行 Scaling 对比，独立生成 `scaling.pdf`，绕过当前 CatMAP `ScalingAnalysis` 的空参数错误。
3. 新增 `descriptor_C_O_H.csv` 和 `plot_C_O_H_volcano.py`。模板每一组催化剂都须填写 **独立、经验证的** `E_C`、`E_O`、`E_H` 逐表面数据。填完后 Notebook 末尾一格绘制 `volcano_C_O_H.pdf`（三页：C–O、C–H、O–H）。
4. **物理定义警示**：当前 `alkene_epoxidation.mkm` 的实际 CatMAP 描述符仍是 `R_s, P_s`；新的 C/O/H 页将当前 R/P 网格中的环氧产物生成速率在每个催化剂的 R/P 坐标处进行数值插值，之后投影到 C/O/H 坐标。故它是*已有模型结果的描述符相关性图*，不是以 C/O/H 为模型输入重新求解后的真实三描述符微观动力学火山图。若需后者，应另构建 R/P/RP 对 C/O/H 的可验证标度关系、输入能量表并重新扫描。
5. `C=-9.28`、`O=-4.37`、`H=-1.11` 为元素原子**固定参考能**，不能直接作为跨催化剂变化的活性描述符。原始 Excel 现有的 OOH、C6H12、O、Ha/Hb 吸附态也不能未经定义就被冒充为 C*、O*、H*。
6. 这次仅进行了 GitHub 代码更新和静态查验；Windows/CatMAP 环境的完整新一轮动态执行仍需在用户本机核实。

在新模型目录运行 `python plot_C_O_H_volcano.py`（需 CSV 15 行全部真实填好），或执行 Notebook 最后一格。若 CSV 未填好，Notebook 会提示缺失的催化剂而非生成虚假火山图。
