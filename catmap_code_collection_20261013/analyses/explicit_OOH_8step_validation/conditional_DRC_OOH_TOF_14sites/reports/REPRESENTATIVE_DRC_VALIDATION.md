# 代表位点 DRC 验证

Ti–Co、Ti–Re、Ti–Mo、Ti–Ni：同一原始八步 CatMAP 模型、统一零库存约束及 h=0.005/0.010 eV。Ti–Mo 本轮新增 44/44 扰动稳态通过；Ti–Co、Ti–Re、Ti–Ni 复用归档的 constrained DRC 与尝试记录，须以 `DRC_NUMERICAL_GATE.csv` 中逐状态门控为准。

官方 `get_rate_control` 不在每次扰动后重解覆盖度，也不保证维持 I1=I2=0；因此按任务要求复用 CatMAP 官方速率与稳态 RHS，并通过已有的七自由覆盖度精确库存坐标求扰动根。所有输出保持 PROVISIONAL。
