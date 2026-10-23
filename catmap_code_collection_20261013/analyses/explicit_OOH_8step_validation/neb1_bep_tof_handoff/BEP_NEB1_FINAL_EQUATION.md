# NEB1 传统 BEP 与八步 CatMAP TOF 衔接

## 同一路径 BEP（本任务定义）

自变量为同一 NEB1 路径反应自由能 `ΔG_NEB1 = G_FS − G_IS`（eV）；因变量为同路径有符号 `ΔG‡_NEB1,signed = G_TS − G_IS`（eV）。使用含截距的一元普通最小二乘，原始负 TS−IS 值不截断：

`ΔG‡_NEB1,signed = 0.354413464088 × ΔG_NEB1 + (0.896960496454 eV)`

配对数据为14个唯一位点。样本内 R²=0.192927，RMSE=0.406917 eV，MAE=0.291653 eV，最大绝对误差=1.202384 eV；逐位点 LOOCV RMSE=0.674522 eV，LOOCV R²=-1.217653。逐位点原始能量、拟合/LOOCV预测及残差见 `BEP_TOF_DATA_MAPPING.csv`。

## 已归档“传统 BEP”的可复用边界

已有包 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization` 的传统 BEP 不是上述严格同路径反应能定义：它用 `EpoxideHHO_strict − Rpair_strict` 作自变量，响应为 `TS1_strict − Rpair_strict`，且排除 Ti–Zr（N=13）。其原归档 α=0.367196621635、β=0.907242240632 eV，R²=0.205428、RMSE=0.417695 eV、LOOCV RMSE=0.713431 eV、LOOCV R²=-1.318030。该模型的自变量定义与样本集合均不同，不能把旧 α/β 用作此处同路径 BEP 系数。旧数据、图件和系数未改动；本次仅对已存14组 IS/TS/FS 配对执行一次符合指定定义的一元拟合。

原始输入：NEB1 配对能量 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_batch_8step\MULTISITE_REACTION_BARRIERS.csv`；旧拟合定义、系数和指标分别见 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization\scripts\run_bep_energy_correction.py`、`H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization\data\BEP_MODEL_COEFFICIENTS.csv`、`H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization\data\BEP_MODEL_COMPARISON.csv`、`H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization\data\BEP_LOOCV_SUMMARY.csv`。

## TOF 对接与解释边界

14位点按 `surface_name` 联结。除 Ti–W 外，TOF 与数值状态沿用 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\conditional_DRC_OOH_TOF_14sites\data\14_SITE_TOF_COVERAGE.csv` / `H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\tiw_tita_targeted_solver_recovery\data\FINAL_14_SITE_STATUS.csv`；Ti–W 使用更新的单点根重评估 `H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\stepwise_kinetics_oxygen_transfer\data\TIW_SINGLE_ROOT_REEVALUATION.csv`，不使用旧汇总表中的失败值，也未重新求稳态。Ti–W 的 log10 值由该校正文件中的精确 TOF 计算，因为校正 CSV 未单独存储对数列。导出同时列出真实 NEB 有符号能量差、拟合值、CatMAP 端点方案有效势垒和八步模型条件性 TOF，各量互不替代。

本式只是这14个 NEB 能量点间的线性样本关系：拟合 R²及误差评价的是 NEB 能量重建，不是 TOF 预测精度。现有 TOF 属八步 CatMAP 数值条件性结果，科学状态仍为 `PROVISIONAL`，且能量映射物理门控含未验证项。因此可描述样本内能垒关联并提出待验证假设，但不能宣称 BEP 已预测催化活性、验证真实位点排序，或将 BEP 预测能垒代入 CatMAP。Ti–W/Ti–Ta 的有符号 TS−IS 为负，均原样保留；拟合不得被解读为这些位点真实活化自由能的物理验证。
