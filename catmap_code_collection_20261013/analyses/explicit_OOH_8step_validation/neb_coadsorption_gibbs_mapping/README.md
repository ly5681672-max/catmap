# NEB共吸附态—独立态 Gibbs 映射核查

本目录的 `rebuild_mapping_audit.py` 是最小只读复现脚本。它读取仓库内已冻结的 DFT_TS_relative、CatMAP energies、现有 Excel 单元格审计及反应能/势垒审计表，用 Decimal 重算端点反应能、差值和 TS 前向保持映射后的逆向势垒。脚本不修改源文件，不调用 CatMAP 稳态求解器，也不运行 DFT/NEB。

在仓库根目录的 Python 环境运行：

```powershell
G:\anaconda\python.exe .\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\neb_coadsorption_gibbs_mapping\rebuild_mapping_audit.py `
  --project-root H:\catmap\catmap `
  --output-dir H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\neb_coadsorption_gibbs_mapping
```

本次运行的逐项输入 SHA256、模型/能量冻结哈希和范围记录见 `SOURCE_RECORDS_AND_HASHES.json`。结论中的数值映射状态不代表结构因果归属；结构同一性、slab/超胞同一性和共吸附能归因仍须独立证据。
