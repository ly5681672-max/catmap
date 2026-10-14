# CatMAP 0.3.1 官方BEP Scaler隔离核验

- 本机版本：0.3.1。
- 使用官方`GeneralizedLinearScaler`的BEP转换矩阵路径；调用`get_coefficient_matrix()`以让官方Scaler按正常调用链初始化约束和矩阵状态。
- 只参数化TS标度关系；元数据确认`ReactionModel.run()`未调用，未生成稳态或TOF。
- TS关系采用初态基准：`E_TS − E_IS = α(E_FS − E_IS) + β`。本项目拟合使用相同有符号定义。

| 路径 | CatMAP α | NumPy OLS α | CatMAP β (eV) | NumPy OLS β (eV) | 样本内 R² | RMSE (eV) | LOOCV RMSE (eV) |
|---|---:|---:|---:|---:|---:|---:|---:|
| NEB1 / path 2 | 0.36719662 | 0.36719662 | 0.90724224 | 0.90724224 | 0.2054 | 0.4177 | 0.7134 |
| NEB2 / path 4 | 0.07921318 | 0.07921318 | 0.19951166 | 0.19951166 | 0.1821 | 0.1172 | 0.1388 |

CatMAP和独立OLS系数在浮点精度内一致。隔离模型配置和调用标志见JSON；模型设置与真实严格运行条件不相混，且没有执行动力学。
