# BEP energy correction optimization

本包以既有13位点CatMAP严格模型能量与NEB1有符号能量差为输入，比较传统BEP、OH/Ha/OOH修正、OH/Ha/OHHa修正和OH调制变斜率BEP。历史单吸附态能量参考尚未完成统一闭合；结果用于模型诊断和探索，不作为新CatMAP基线。

## 目录

- `data/`：只读输入快照及生成的模型、系数、误差、消融、共线性与TS1重建表。
- `scripts/`：可复现分析入口。
- `figures/`：七类图件，PNG/PDF各一份。
- `origin_exports/`：逐位点数据、拟合/残差、系数、曲线及消融数据。
- `reports/`：中文分析报告、SHA256清单及运行日志。

## 复现

从本目录运行：

```powershell
python .\scripts\run_bep_energy_correction.py
```

脚本仅读取本包 `data/` 快照并写入本包；不调用CatMAP稳态求解，不改写任何源项目文件。

11位点诊断集排除Ti-W与Ti-Ta。报告分别列出11位点重新拟合以及13位点系数应用到同一11位点的结果，两者不可混同。
