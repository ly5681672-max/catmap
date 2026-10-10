# lyq_alkene_epoxidation1：原版八步框架 + BEP 描述符筛选 + CatMAP 活性火山图

**本项目与 `lyq_alkene_epoxidation` 同构；反应网络仍然是原来的八个基元步骤。**
本次修复纠正了此前误将新版换成 R/P/RP 三步共吸附模型的偏差。
原项目不变，所有 DFT 原始总能及 Excel 不变。

## 1. 两类数据不可混淆

- `1-generating_input_file/data_only_not_executable.py`：已从原版**逐字恢复**完整 Ti 原始总能、`ref_dict`、形成能函数、频率字典和 CatMAP 表格解析逻辑。该脚本名虽写 `not_executable`，但其原本含有实际输出语句；不要在不审阅时重复执行。
- `1-generating_input_file/generate_input1.py` 以及每个 `ti*_energies.txt/.csv`：已恢复原版代码与原版逐催化剂输入。
- `2-creating_microkinetic_model/energies.txt` 与 `energies_origin.txt`：已恢复原版**8步机理所需的完整形成能**。不能改成只含 R/P/RP 的数据表。
- `1-generating_input_file/energy_audit.csv`：保留之前直接从 Excel 提取的完整 NEB 反应 IS/marker/FS、相对能及暂定有效势垒；作为**独立 BEP 筛选的数据源**，不在 CatMAP 中冒充第4步的 TS。
- `己烯环氧化.xlsx`：原始工作簿保留不变。

`ref_dict['C']=-9.28`、`ref_dict['O']=-4.37`、`ref_dict['H']=-1.11` 为固定**原子参考能**，不是随催化剂变化的描述符。

## 2. C/O/H-family 的实际候选描述符

- `C6H12_ads` = \(G_f(C_6H_{12}^*)-G_f(C_6H_{12,g})\)；代表 C 相关烯烃吸附。
- `O_form` = \(G_f(O^*)\)；代表 O*。
- `Ha_form`、`Hb_form` = \(G_f(Ha^*)\)、\(G_f(Hb^*)\)；分别代表两种 H 吸附位置。
- `OH_form` 等可作为**扩展**描述符，但默认筛选严格限定 C/O/H-family。

环氧化的目标是来自 Excel 的**完整净反应**暂定正向有效势垒 \(\max(0,\Delta G_{marker},\Delta G_{FS})\)；这不是 NEB 全部图像经过 saddle 验证的真实势垒，也不是原八步机理第4步的精确势垒。
因此回归应称为 **BEP-like descriptor–barrier linear relation**，不能直接将其宣称为已验证的传统同一步 BEP 定律。

默认以原版 .mkm 采用的**14 个催化剂**做筛选（原版排除 TiNb）；分析文件中保留 TiNb 原始数据，未删除。用线性 OLS + 留一交叉验证（LOOCV）筛选两描述符，避免只按训练 R² 排序。

本轮筛选的 C/O/H-family 最优组合为 **O_form + Ha_form**：14组 R²≈0.454、LOOCV RMSE≈0.156 eV。
若允许扩展描述符，**OH_form + Ha_form** 更强：14组 R²≈0.679、LOOCV RMSE≈0.114 eV。因用户当前指定 C/O/H-family，.mkm 仅采用 **O_s / Ha_s** 两坐标。

## 3. CatMAP 火山图与 BEP 图的不同含义

`2-creating_microkinetic_model/alkene_epoxidation.mkm` **8步反应、气相活度、温度、物种定义均按原版**；仅改为：
```python
descriptor_names = ['O_s', 'Ha_s']
descriptor_ranges = [[-4.0, 2.5], [-4.0, 0.8]]
scaling_constraint_dict = {
    'O_s': ['+', 0, None],
    'Ha_s': [0, '+', None],
}
```

`1-generating_input_file/bep_descriptor_analysis.py` 输出：
- `bep_screening.csv`、`bep_predictions.csv`：候选线性拟合、LOOCV 与逐表面残差；
- `bep_fit.pdf`：实际暂定势垒 vs 拟合势垒；
- `bep_barrier_map.pdf`：拟合的二维能垒平面，**不是 TOF 火山图**。

`2-creating_microkinetic_model/test.ipynb` 将先运行上述筛选，然后运行**原版八步 CatMAP**，输出标准 rate / all_rates / coverage / scaling 图，并额外输出：
- `production_rate_volcano.pdf`：真正基于 CatMAP `production_rate_map` 的 O*/Ha* 二维环氧产物生成速率图，带催化剂位置标注；
- `production_rate_table.txt`：完整数值表。

**模型边界**：在八步网络中，CatMAP 对其他吸附态和 TS 使用 generalized linear scaling，自动拟合的 TS 能量不会严格等于 BEP 筛选值；且原始第4步使用的名义 TS 与整体共吸附 NEB 的化学阶段不完全一致。因此火山图只能作为原机理假设下的**探索性活性映射**，不能直接拿来作为验证的 TOF 预测。要得到严格 BEP 控制的第4步，必须额外提供该基元步骤对应的真实 NEB 路径并建立能量一致的过渡态约束。

## 4. 本地运行

```powershell
cd H:\catmap\catmap
git status --short
git pull --ff-only origin main
```

在 `lyq_alkene_epoxidation1/2-creating_microkinetic_model/test.ipynb` 中重启 Kernel，依次运行所有单元。第一单元拟合并生成 BEP 图，第二单元求解八步 CatMAP，第三单元输出原始类的速率/覆盖度/Scaling PDF，第四单元生成 `production_rate_volcano.pdf`。

单独检验描述符：
```powershell
python lyq_alkene_epoxidation1/1-generating_input_file/bep_descriptor_analysis.py
```

该脚本仅写本项目 2-creating_microkinetic_model 文件夹的 BEP 报表/诊断图，不修改任何 `energies.txt`、Excel、`.mkm`。

目前已完成 GitHub 源码级核查及 14/15 组候选回归数值核验。更新后的八步 CatMAP 全网格在用户 Windows 环境是否全部收敛、全部图片是否生成，仍需运行时检查。
