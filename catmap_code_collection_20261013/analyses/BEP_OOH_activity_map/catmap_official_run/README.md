# 四参数 OOH* 扩展 BEP：CatMAP 0.3.1 正式稳态复核

本子目录用于隔离复核四参数扩展 BEP 对 TS1 的替换影响，并通过本机 CatMAP 0.3.1 `ReactionModel` / `ThermodynamicScaler` / `SteadyStateSolver` 计算三组13位点TOF和 Ti–Co 条件参照二维网格。上级目录中的历史解析活性图和数据保持不变。

## 目录

- `scripts/run_bep_ooh_activity_map.py`：生成输入快照、复拟合精确BEP、执行原始基线/全样本/LOOCV及5×5、15×15、61×61真实CatMAP节点、验收并出图。
- `inputs/snapshots/`：只读来源输入快照与哈希/描述符范围来源。
- `inputs/cases/`：自动生成的每次独立CatMAP输入；大批量文件不纳入Git。
- `runs/`：CatMAP节点的结果、日志、断点缓存；不纳入Git。
- `data/`：可审核汇总CSV、系数、覆盖度、六步净速率、节点状态和SHA256清单。
- `figures/`：TOF对比、TS1重建误差和CatMAP二维活动图，PNG/PDF/SVG。
- `reports/`：阶段日志与最终科学报告。

## 复现

在PowerShell运行（环境变量只影响该进程）：

```powershell
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
G:\anaconda\python.exe -X utf8 .\analyses\BEP_OOH_activity_map\catmap_official_run\scripts\run_bep_ooh_activity_map.py --stage all --python G:\anaconda\python.exe
```

脚本通过SHA256复用输入完全相同且已验证的缓存；`--stage baseline`、`--stage fit`、`--stage grid --resolution 5|15|61`可分阶段断点续算。最终61×61表逐个记录 `CONVERGED`、`CONVERGED_ADAPTIVE`、`REUSED_CACHE` 或 `FAILED`；失败不填零、不插值。

## 物理与数值口径

- `OOH*` 是历史单吸附态 formation energy，作为BEP外部描述符；不加入六步机理，也不声称与桥连OOH共吸附态同源。
- BEP输出是有符号 `G(TS1)-G(Rpair)`。预测TS低于初态/终态时，保留预测值并标记；CatMAP 0.3.1速率模板按三态能量最大值计算有效正逆势垒（端点控制），不在输入层人为修正TS。
- 所有模型仍为333.15 K、原六步机理、原严格有效储库浓度、`tolerance=1e-70`、75位mpmath及单有效中心粗粒化定义。
- 二维图是Ti–Co固定能量背景下的条件活性图；叠加的13个位点符号是其真实描述符坐标，不代表13个点均在图中背景条件下重新计算。
