# 实际运行记录

- 输入：`H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021\`
- Python：`G:\anaconda\python.exe`
- CatMAP：0.3.1（通过项目治理预检；正式八步模型解析配置保持原样）
- 执行：对归档的11个已接受覆盖度根调用 CatMAP parse/compile、官方稳态函数、速率常数、净反应速率及 TOF；无 `ReactionModel.run()`、无求根。
- 结果：11个根均成功完成只读重新评价；另3个位点无 accepted-root JSON，标为无可复核根。
- 生成文件：`saved_root_local_flux_audit.csv`、`data/THREE_SITE_REEVALUATION.csv`、`data/OLD_VS_NEW_QUALITY_STATUS.csv`。
- 代码校验：`G:\anaconda\python.exe -m py_compile evaluate_saved_eight_step_roots.py numerical_gate_audit_20261021/scripts/build_gate_tables.py` 通过；两脚本实际运行成功。
- 原始 `.mkm`、能量文件及历史结果没有修改；无 Git commit/merge/push。

