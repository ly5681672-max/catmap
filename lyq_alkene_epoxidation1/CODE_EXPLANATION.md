# NEB 共吸附独立对照版说明

与原目录的子目录、文件名称和相对位置保持相同。**此副本**对应整体共吸附反应；旧项目仍为八步基元机理。二者的 `energies.txt` 不能混用。

## 运行

1. 激活可导入 `catmap` 的 Python 环境。
2. `cd lyq_alkene_epoxidation1/2-creating_microkinetic_model`。
3. `python -c "from catmap import ReactionModel; ReactionModel(setup_file='alkene_epoxidation.mkm').run()"`。
4. 在 `test.ipynb` 中复现和检查。

## 生成数据

在 `1-generating_input_file` 执行 `python generate_input1.py`；其输入直接来自本目录上一层的 Excel，写入同目录的 `energy_audit.csv` 及模型目录的 `energies.txt`。额外依赖 `artifact_tool`，如未安装可直接使用已提交的能量表。

## 反应与物种

- `R*`: (H2O2)(C6H12) 完整共吸附态，C6H14O2。
- `P*`: (C6H12O)(H2O) 完整共吸附态，C6H14O2。
- `RP*`: 按 `max(G_IS,G_marker,G_FS)` 建立的**数值有效 TS**；未验证为 CI-NEB 最高点。
- `* + H2O2_g + C6H12_g <-> R*`
- `R* <-> RP* -> P*`
- `P* <-> C6H12O_g + H2O_g + *`

气相后缀的 pressure 参数在此仅为液相近似有效活度；TOF 尚未通过 CatMAP 实跑及完整 NEB 校准。不同自由能修正基准不可混用。TiMo 的末态取 Excel 第21行组合总能，TiTa 的 G 列单独偏移。

## 与原版直接比对

- 原：`lyq_alkene_epoxidation/2-creating_microkinetic_model/alkene_epoxidation.mkm`
- 新：`lyq_alkene_epoxidation1/2-creating_microkinetic_model/alkene_epoxidation.mkm`

项目根目录及子目录命名仅差一个 `1`，其他主要对应文件路径保持一致。
