# CatMAP 项目快速导览：按微观动力学分析流程阅读

归档位置：`H:\catmap\catmap\catmap_code_collection_20261013`  
本指南盘点日期：2026-10-14  
用途：帮助快速找到每个研究阶段的输入、模型、计算结果、验证证据和图表；它是导航文件，不替代各分析包的原始报告。

## 先看结论：从哪里开始

如果只花十分钟熟悉项目，建议按下面顺序打开：

1. `PROJECT_ARCHITECTURE.md`：归档规则、原始数据保护与目录约定。
2. `catmap_official_run/alkene_epoxidation.mkm` 和同目录 `energies.txt`：了解六步机理、物种组成和CatMAP TableParser输入。**注意根目录模型是静态验证配置，不是稳态求解入口。**
3. `catmap_official_run/results/official_strict_tolerance/`：逐催化位点的严格容差官方模型副本；主要基线结果汇总见第4节所列机制数据包。
4. `analyses/microkinetic_mechanism_13sites/reports/MECHANISM_ANALYSIS_REPORT.md`：查看稳态覆盖度、通量、有效势垒、DRC、数值边界和13个位点汇总。
5. 根据问题进入 `analyses/` 下的专题包；其 `reports/` 通常先给结论和限制，`data/` 给逐点数据，`scripts/` 给生成流程，`figures/` 给图，`origin_exports/` 给Origin表格。
6. 查找单个文件时使用同目录的 `CATMAP_FILE_CATALOG.csv`；它按相对路径列出归档根目录内的全部文件、类型、大小和用途线索。

## 研究主线：从原始能量到机理解释

```text
原始Excel/NEB与Gibbs记录（只读来源）
        ↓ 状态、催化剂、单元格与能量参考审计
CatMAP energies.txt + 六步 .mkm（统一输入层）
        ↓ CatMAP 0.3.1 官方 Parser / Scaler / Solver
严格稳态：TOF、覆盖度、反应净通量、有效势垒、DRC
        ↓ 独立复核与按 surface_name 匹配
机制分析 → 描述符/标度关系 → 官方描述符网格或经验关系
        ↓
图件、Origin数据、报告、哈希/测试证据
```

建议分析任何新问题时，先回答四件事：

- **输入是什么？** 记录使用哪份 `energies.txt`、`.mkm`、温度、储库配置和CatMAP版本。
- **计算属于哪一层？** 静态Parser/反应能检查、官方稳态求解、官方描述符网格，还是只对既有数据作统计拟合。
- **数据如何连接？** 跨文件按 `surface_name` 等明确键一对一连接，不按表格行序拼接。
- **结论证据等级是什么？** 区分原始DFT/NEB、官方CatMAP结果、Scaler生成的虚拟点/VectorMap插值、经验回归和条件性参考转换。

## 1. 根目录与通用目录

| 路径 | 内容与用途 | 阅读/使用提醒 |
|---|---|---|
| `PROJECT_ARCHITECTURE.md` | 归档结构、只读数据原则、复制/校验规范 | 新成果按此结构归档；不要散放根目录或覆盖历史包。 |
| `README.md` | 汇总来源、复制范围、文件类型统计和归档限制 | 这是2026-10-13的汇总说明；若后续归档新增内容，统计数字不会自动更新。 |
| `catmap_official_run/` | 官方CatMAP输入、逐位点模型、运行结果和官方图 | 根目录 `.mkm` / `run_catmap.py` 是静态验证入口；稳态入口另在逐位点目录。 |
| `models/v1_serial/` | V1串行NEB模型的映射、物种和能量验证脚本/配置 | 属于模型阶段档案，不能不核版本地替代当前六步严格模型。 |
| `scripts/` | 项目级Excel/NEB提取、CatMAP版本审计、阶段报告构建脚本 | 读取脚本中的绝对路径；归档不等于已路径迁移或可移植运行。 |
| `tests/` | Excel提取、计量守恒、CatMAP能量语义、串行网络等回归测试 | 测试只对其自身夹具/代码范围有效，不自动证明所有物理假设成立。 |
| `v2_global_reference/` | V2全局Gibbs参考、slab、NEB源行映射与热力学网络审计代码 | 独立审计支线，不是严格稳态基线；`data/construct_global_gibbs.pre_source_row_fix.py` 明确是修复前快照，不可当权威映射。运行前先查本支线报告/状态与输入版本。 |
| `presentation_revision_20261013/` | 2026-10-13汇报PPT图示修改脚本 | 是演示文稿处理工具，不参与化学能量或CatMAP计算。 |
| `analyses/<专题名>/` | 独立派生分析包，含数据、脚本、图、Origin导出和报告 | 每个研究问题独立保存；细目见下文及文件索引。 |

## 2. CatMAP 0.3.1模型与计算入口

### 标准输入

| 文件 | 作用 | 阅读重点 |
|---|---|---|
| `catmap_official_run/energies.txt` | CatMAP `TableParser`表格：surface、site、species、formation energy、频率、reference | 这是六步网络的形成能输入快照。不要把旧单物种表或局部NEB能差直接混入此表。 |
| `catmap_official_run/alkene_epoxidation.mkm` | 六步反应、物种composition/n_sites、温度/描述符、储库与Scaler/Solver配置 | 文件头写明 `Conditional static-validation setup`，`surface_names=['tire']`，`resolution=1`，脚本不会调用稳态运行。它适合解析、守恒和Ti–Re静态势垒核对，不是14位点正式运行配置。 |
| `catmap_official_run/run_catmap.py` | 对14个surface逐一调用官方 `parse()`/`verify()`；Ti–Re另核反应能和端点势垒 | 脚本明确 `ReactionModel.run()` 未调用；不能把它的PASS说成已求稳态。 |

六步机理为：反应物共吸附 → NEB1显式TS1 → 环氧产物脱附 → NEB2显式TS2 → 表面水形成 → 水脱附。第2/4步保留显式TS1/TS2；第1/3/5/6步在探索性动力学模型里使用CatMAP端点势垒近似。Rpair/TS1/EpoxideHHO以C6H14O2 composition表示；HHO/TS2/OHHa/H2O以H2O总组成表示。`n_sites=1`是单有效中心粗粒化约定，不代表真实几何占位数已单独证明。

### 逐位点实际模型/结果层

| 路径 | 用途 |
|---|---|
| `catmap_official_run/results/official_strict_tolerance/<surface>/` | 逐位点严格容差模型文件副本，含 `.mkm`、`energies.txt`、`run_catmap_official.py`，部分位点另有自由能图。是识别严格基线运行配置的首要位置。 |
| `catmap_official_run/results/official_14/<surface>/` | 第一轮14位点官方模型运行档案；读结果时应与严格容差版本区分。 |
| `catmap_official_run/results/official_14_retry_damping50/` | 初轮失败位点的阻尼参数重试快照；属于数值重试，不应覆盖严格基线。 |
| `catmap_official_run/results/tire_official/` | Ti–Re官方运行的单独快照，需结合内部日志/config判定其是静态还是动力学运行。 |
| `catmap_official_run/results/activity_analysis_round1/` | 复用既有稳态结果的活性、覆盖度、ΔG3与DRC分析及图。 |
| `catmap_official_run/results/official_plots_round2` 至 `round5` | CatMAP官方绘图的多轮输出；round4中含先前图的归档副本，不要重复计为新结果。 |
| `catmap_official_run/TS1_TS2_descriptor_grid/`、`catmap_official_run/O_H_to_TS1_TS2_to_TOF/` | 官方描述符网格脚本、网格结果缓存/日志或绘图数据。具体哪个是原始网格、哪个是VectorMap插值，需看网格摘要和各专题报告。 |

**重要数值边界：**机制专题报告记录13个位点严格基线均收敛，Ti–Zr仍未验证/不纳入13点分析；温度333.15 K。第1/3/5/6步的势垒是端点近似而非独立DFT过渡态。方向速率接口在近零空位处有双精度局限；遇到方向速率异常，应优先引用高精度净通量、覆盖度和TOF的独立检查，不把负方向输出直接当真实单向通量。详细依据见 `analyses/microkinetic_mechanism_13sites/reports/MECHANISM_ANALYSIS_REPORT.md`。

## 3. 分析专题包导航

以下七个包均以“一个科学问题一个目录”组织。括号中是盘点时各包文件数（可能随后续工作变化）；准确的逐文件路径见 `CATMAP_FILE_CATALOG.csv`。

### A. `analyses/microkinetic_mechanism_13sites/`（132个文件）

**问题：** 13个已收敛位点的稳态机制如何不同？  
**先读：** `reports/MECHANISM_ANALYSIS_REPORT.md`。  
**核心数据：** `data/convergence_13sites.csv`、`coverage_13sites.csv`/`coverage_13sites_long.csv`、`reaction_energy_13sites.csv`、`effective_barriers_13sites.csv`、`directional_reaction_rates_13sites.csv`、`TOF_flux_consistency.csv`、`rate_control_13sites.csv`/`rate_control_matrix.csv`、`explicit_TS_signed_energy_differences.csv`。  
**复现脚本：** `scripts/analyze_mechanism_outputs.py`（分析/图表），`scripts/run_missing_official_drc.py`（补充官方DRC计算，运行会触发CatMAP）。  
**输出：** 覆盖度、正逆/净速率、势垒、DRC热图与Origin长表/矩阵。该包目前还有 `data/official_drc_runs_supplemental/<surface>/`，保存补充运行配置、日志、pickle与rate_control CSV。  
**注意：** 报告结论优先于仅看图；DRC是局部敏感性，不等同唯一决速步。方向API的已知精度问题和低TOF警示在报告中有说明。

### B. `analyses/dG3_coverage_DRC_mechanism/`（31个文件）

**问题：** 第3步标准/有效ΔG、覆盖度、DRC与TOF的关联，以及TS1–TS2官方网格峰附近的响应。  
**先读：** `reports/MECHANISM_CORRELATION_REPORT.md`、`reports/TIFE_DATA_AND_MECHANISM_AUDIT.md`、`reports/EXISTING_WORK_REUSE_AUDIT.md`。  
**核心数据：** `data/combined_dG3_coverage_DRC_TOF.csv`、`data/dG3_coverage_TOF_correlations.csv`、`data/TS1_TS2_peak_mechanism_data.csv`、Ti–Fe审计表和输入哈希。  
**复现脚本：** `scripts/build_mechanism_package.py`；图/Origin文件在对应目录。  
**注意：** ΔG3–TOF相关部分由第3步端点势垒定义耦合；Ti–Fe slab来源偏移尚未从结构证据上闭合；虚拟网格峰没有峰区DRC直接证明。

### C. `analyses/TS1_TS2_descriptor_volcano/`（39个文件）

**问题：** TS1_s + TS2_s 描述符平面的经验响应与CatMAP官方网格是否表现内部活性峰。  
**先读：** `reports/TS1_TS2_VOLCANO_REPORT.md`、`ARCHIVE_README.md`、`reports/EMPIRICAL_NOTES.md`。  
**脚本：** `scripts/run_official_grid.py`（官方网格，可能启动CatMAP计算）、`scripts/probe_official_model.py`（官方模型/接口检查）、`scripts/generate_empirical.py`（经验二次响应面）。  
**结果：** `figures/catmap_official/` 是CatMAP网格和VectorMap图；`figures/empirical/` 是统计响应面；`origin_exports/catmap/` 与 `origin_exports/empirical/` 分开保存数据。  
**关键区别：** 官方50×50网格节点是CatMAP计算；400×400 VectorMap是对稀疏网格的官方映射/插值，不是160,000次独立稳态求解。经验二次图不属于CatMAP动力学计算。实际催化剂锚点和标度虚拟点也必须区分。

### D. `analyses/O_H_to_TS1_TS2_to_TOF/`（61个文件）

**问题：** 旧O/Ha/Hb特征能否经验重建TS1/TS2，并连接到TOF与既有CatMAP TS1–TS2网格。  
**先读：** `reports/O_H_to_TS1_TS2_TOF_REPORT.md`、`data/TS1_TS2_VOLCANO_REPORT.md`、根目录 `README.md`。  
**脚本：** `scripts/analyze_oh_to_ts_to_tof.py`、`scripts/write_sha256_manifest.py`。  
**数据：** `data/matched_data.csv`按surface_name匹配；`network_state_energy_long.csv`、TS预测/BEP/回归指标、`reference_slab_comparison.csv`分别记不同分析层。`catmap_official_run/`是官方网格证据快照。  
**图件与导出：** `figures/scaling/`、`figures/bep/`、`figures/tof_maps/`；`origin_exports/maps/`包含经验模型网格和已保存CatMAP VectorMap插值网格。  
**注意：** O/Ha/Hb为历史单物种参考，和六步模型的共吸附态不是同一个定义；slab及热修正一致性尚有限制。包内同时存在经验响应面与CatMAP已有插值，两者文件名已区分，引用时勿混为一类。

### E. `analyses/empirical_O_H_descriptor_analysis/`（34个文件）

**问题：** 历史O*/Ha*/Hb*形成能与严格CatMAP真实TOF的经验关系。  
**先读：** `reports/NOTES.md`、`reports/FULLGRID_NOTES.md`。  
**脚本：** `scripts/analyze_legacy_descriptors.py`、`scripts/export_fullgrid_for_origin.py`。  
**数据：** `data/matched_data.csv`、`fit_metrics.json`、`reference_slab_comparison.csv`、`workbook_source_audit.csv`。  
**图与Origin：** O/Ha/Hb单变量经验图；O+Ha、O+Hb经验二维图及其全矩形200×200、400×400网格、催化剂点和峰值表。  
**注意：** 属经验拟合，不是CatMAP连续网格；旧单物种能量相对当前严格网络的物理slab参考一致性未完全证明。

### F. `analyses/C_O_H_to_TS1_scaling/`（29个文件）

**问题：** 核查独立C*是否存在，并比较可用O/Ha/Hb特征对TS1的线性标度。  
**先读：** `reports/CSTAR_SOURCE_AUDIT.md`、`reports/TS1_DESCRIPTOR_SCALING_REPORT.md`、`reports/VALIDATION_LOG.md`。  
**脚本：** `scripts/analyze_ts1_descriptors.py`。  
**数据：** 来源扫描/工作簿审计、匹配数据、拟合指标/逐催化剂残差、C*模板与哈希。  
**结论边界：** 该包审计中没有找到真实独立C*，没有用C6H12*冒充C*；O/H样本内拟合较弱且历史参考有限。后续含碳反应物/中间体探索见G包。

### G. `analyses/Carbon_containing_descriptor_to_TS1/`（57个文件）

**问题：** 用已有含碳吸附物种与O/H数据探索TS1形成能和NEB1有符号能量差的重建。  
**先读：** `reports/ANALYSIS_REPORT.md`。  
**脚本：** `scripts/analyze_carbon_ts1_scaling.py`。  
**候选数据：** `data/descriptor_source_inventory.csv`、`data/descriptor_availability.csv`、全组合指标、系数和逐点残差、workbook crosswalk。  
**图/Origin：** TS1与NEB1拟合对照、含碳/O/H性能图、完整矩形经验响应面和200×200/400×400 XYZ/矩阵。  
**注意：** 历史C6H12*、C6H12O*、Rpair等真实物种不能混称C*；某些含产物/反应前状态的组合具有事后或网络耦合性质。历史参考与当前严格模型参考之间仍有条件性限制；图是经验回归，不是官方CatMAP网格。

## 4. 如何读机制数据，而不把不同量混在一起

| 数据 | 物理/计算含义 | 常见误读 |
|---|---|---|
| `formation_energy` | 相对输入参考态/原子参考的物种能量 | 不等于吸附能，跨不同slab/工作簿需先核参考。 |
| 标准 `ΔG` | 仅由CatMAP输入能量和反应化学计量得到的反应能 | 不含指定储库化学势影响。 |
| 有效 `ΔG` | 指定温度、储库与模型配置下的条件反应能 | 不能与标准ΔG当作同一量。 |
| 显式NEB `TS−IS` | 输入NEB态相对其初态的有符号差值 | 原始负差保留；不要事后截成0，也不要把端点反应能叫活化能。 |
| 有效正/逆势垒 | CatMAP按实际反应式、能量、储库/端点模式给出的模型势垒 | 第1/3/5/6步是端点近似，不是独立DFT TS。 |
| 覆盖度 | CatMAP稳态占位分布 | 最大覆盖物种不自动等于唯一速率控制步骤。 |
| DRC | 对指定状态能量的局部速率敏感度 | 不是全局唯一控速步骤；极低TOF或双精度方向速率异常要审慎。 |
| 官方二维网格 | CatMAP在描述符坐标下调用求解器形成的虚拟样本 | 网格虚拟催化剂不是新DFT催化剂；VectorMap高分辨率图可能只是对稀疏网格的插值。 |
| 经验拟合/响应面 | 对已有催化剂数据做统计回归并在指定矩形评估 | 不得标注为CatMAP连续微观动力学结果；矩形外/无真实点区域是模型延拓。 |

## 5. 按研究问题选择最短路径

- **“某位点为什么活性不同？”** 先用 `microkinetic_mechanism_13sites` 的覆盖度、六步净通量、有效势垒、DRC和报告；需要看ΔG3耦合时再读 `dG3_coverage_DRC_mechanism`。
- **“TS1/TS2是否能由简单能量描述？”** 先看 `Carbon_containing_descriptor_to_TS1`（含碳候选扩展），再看 `C_O_H_to_TS1_scaling` 的独立C*来源审计及 `O_H_to_TS1_TS2_to_TOF`。
- **“描述符空间的CatMAP活性峰在哪？”** 先看 `TS1_TS2_descriptor_volcano` 官方网格报告及数据，再区分经验二次面与官方网格；确认锚点、Scaler约束、TS标度模式和插值层级。
- **“图中坐标/峰值/误差从哪来？”** 找专题包 `origin_exports/` 与 `data/`，再沿脚本回溯；不要从图像像素读数代替CSV。
- **“如何复现？”** 先检查脚本中的输入绝对路径、CatMAP版本、`results_directory`/输出位置和缓存策略；复现前复制到隔离目录，禁止直接覆盖历史缓存或正式基线。
- **“遇到能量参考冲突？”** 从当前 `.mkm` 实际指定的 `energies.txt` 追溯来源单元格、clean slab、元素参考和热修正；V2支线脚本仅作为静态审计工具，不能默认为权威输入转换。

## 6. 归档、复现与新增工作的约定

1. 原始Excel、DFT/NEB文件、正式CatMAP输入和已验收基线只读；新分析建新专题目录，不覆盖旧图/CSV/cache。
2. 新专题建议结构：`analyses/<topic>/{data,figures,origin_exports,scripts,reports}`。官方模型/求解副本放 `catmap_official_run/`，通用脚本和测试分别放根目录 `scripts/`、`tests/`。
3. 原始/派生输入用清晰来源名和SHA256记录。按 `surface_name`、物种名、NEB路径/角色等科学键连接；输出保留来源，不以行号连接。
4. 归档脚本可能保留C盘绝对路径。复制到H盘不等于已改写路径，也不证明可独立运行。重新运行前先审阅脚本与配置，选择隔离输出目录。
5. 对CatMAP官方图，记下本机版本、scaler、solver、约束、descriptor ranges/resolution、是否真实求解或VectorMap插值。对经验图，明确拟合模型和有效数据范围。
6. 每个新分析包至少保存输入清单、结果报告、运行日志/测试结果和哈希清单；失败结果也保留，不能用“TOF=0”代替未收敛。

## 7. 归档现状快照与证据限制

本次清点的根目录主要组成：`analyses/` 7个专题包；`catmap_official_run/`含标准输入、运行结果、绘图轮次和TS1/TS2网格；另有 `models/`、`scripts/`、`tests/`、`v2_global_reference/`、`presentation_revision_20261013/`。各目录和文件的全量相对路径、字节数、扩展名、用途说明见 `CATMAP_FILE_CATALOG.csv`。

此指南依据现存归档的README、架构约定、模型配置和专题报告归纳。它不表示在H盘独立重放过所有脚本，也不表示归档内所有历史实验结果都属于当前正式基线。若报告与缓存/图件出现版本差异，以来源可追溯的输入和报告中的阶段标记为准；不清楚时先做只读审计。
