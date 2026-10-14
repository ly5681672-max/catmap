# BEP Scaling and TS Reconstruction

本专题复核CatMAP六步己烯环氧化模型的NEB1/NEB2 BEP关系、额外共吸附/中间体能量项、过渡态能量重建及其与已保存CatMAP TS1–TS2 VectorMap的对应。此包是探索性分析和归档材料，不是新的正式CatMAP基线。

## 目录

- `data/`：从H盘归档只读复制的原始输入快照、官方Scaler系数与分析摘要。
- `scripts/`：CatMAP官方BEP静态核验、OLS/LOOCV表格生成、绘图、报告和校验脚本。
- `figures/`：PNG/PDF科研图。
- `origin_exports/`：完整逐催化剂数值、BEP拟合、残差、TS重建和VectorMap近邻映射CSV。
- `reports/`：中文结果、数据审计、Scaler核验、源输入哈希、校验日志和文件SHA256。

## 复现

在Windows PowerShell运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\run_bep_package.ps1
```

脚本使用本机CatMAP 0.3.1、NumPy/Pandas/scikit-learn/Matplotlib。官方BEP probe仅构建系数矩阵；严禁调用`ReactionModel.run()`。绘图背景只读取既有的400×400 CatMAP VectorMap插值结果，不进行新网格求解。

## 科学边界

NEB正向/逆向原始有符号能差完整保留；负的TS−IS差值不截断。历史单吸附态能量与当前严格共吸附网络态的slab/热修正参考尚未闭合，扩展项只作探索性诊断。`ΔG1`是气相反应物到Rpair的整体自由能差，不等于单独定义的共吸附相互作用能。BEP/TS重建不能替代新催化剂的真实NEB或正式稳态计算。
