# 验证日志

- 运行入口：`G:\anaconda\python.exe scripts\analyze_ts1_descriptors.py --project-root H:\catmap\catmap\catmap_code_collection_20261013 --raw-project-root H:\catmap\catmap\lyq_alkene_epoxidation --workbook H:\Doctoral_research\己烯环氧化.xlsx --workbook H:\catmap\catmap\lyq_alkene_epoxidation\己烯环氧化.xlsx`
- 输入条件：仅读取现有 `matched_data.csv`、`energies.txt`、来源/参考审计、旧转换表/脚本及两份 Excel；不调用 CatMAP，不计算 TOF，不运行 DFT。
- 关键断言：13 个 `surface_name` 一对一匹配；全部标记 `CONVERGED` 且温度 333.15 K；`tizr` 排除；TS1 键无重复；所有拟合输入有限。
- 源数据核对：两份 Excel 各 27 张工作表的显式 C* 标签数均为 0；旧形成能表独立 C* 物种数为 0；没有计算 C 相关回归。
- 回归验证：全部 3 个可用单变量及 3 个可用 O/H 双变量输出 R²、RMSE、MAE、最大绝对误差、LOOCV RMSE、原单位系数及逐位点残差。
- 复跑检查：以 `descriptor_fit_metrics.csv` 和 `descriptor_fit_predictions.csv` 的 SHA256 比对连续两次运行；成功条件为两者均逐字节一致。
- 绘图检查：性能图及 3 张双描述符实际-预测图均输出 PNG/PDF，所有 13 个真实催化剂名称均由对应点表标注；不生成 C* 网格/响应面。
- 本日志记录执行设计；最终逐文件 SHA256 见 `data/SHA256SUMS.csv`。
