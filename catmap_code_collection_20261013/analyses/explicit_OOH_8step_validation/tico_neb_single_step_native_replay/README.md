# Ti–Co CatMAP 原生 NEB 单步复现

此目录为隔离的单步数值验证，不属于正式八步模型，也不计算稳态、TOF 或 DRC。

两个反应将完整 NEB 的 IS、TS、FS 分别表示为同一化学总组成的三个“状态标签”。模型中的 `n_sites=1` 仅为 CatMAP 单反应解析/位点平衡所需的等位状态账本占位约定，**不是** DFT 几何对真实共吸附态占据一个位点的证明。该小模型只验证路径上的 Gibbs 差、CatMAP 正逆有效势垒和官方速率常数接口。

运行：

```powershell
G:\anaconda\python.exe .\scripts\replay_catmap_neb.py
```

运行时目录应为此目录。输出写入 `results/`。正式八步 `.mkm`、原始 `energies_dft_neb.txt` 和既有14位点结果均不被修改。
