# 显式 OOH 八步模型状态

## 结论

原始八步反应式已由本机 CatMAP 0.3.1 实际加载并解析为 8 个基元反应、9 个稳定吸附物种、2 个过渡态；八步机理中保留稳定 OOH*。原始 Ti–Co Excel 六个 NEB Gibbs 值及其公式缓存、现有 CatMAP 输入中的数值映射可复现，但**NEB 与独立吸附态的共同 slab/结构/热修正等价性仍未证实**。因此当前不能宣布 NEB1/NEB2 已通过物理映射验收，也未达到正式全局 TOF/DRC 分析门槛。

## 官方解析事实

- CatMAP：`0.3.1`，安装位置 `G:\anaconda\Lib\site-packages\catmap`。
- 原始模型：`lyq_alkene_epoxidation/2-creating_microkinetic_model/alkene_epoxidation.mkm`。
- 8 个解析基元反应与 9 个稳定吸附物种见 `../data/eight_step_parse_screen_summary.json`。
- CatMAP 将 `OOH-C6H12_s` 与 `Hb-O_s` 解析为过渡态；它们不进入稳态覆盖度守恒式。
- 9 个稳定吸附物种均由原始 `.mkm` 与 CatMAP Parser 得到 `n_sites=1`。这只是模型记账占位，不是几何位点数的独立实证。

## NEB 映射状态

NEB1 数值对应八步候选第 4 步，NEB2 数值对应第 6 步；两条路径的 G 数值来源已核验，但这不等于共吸附 NEB 的完整局部环境已由现有独立物种表表达。两组 NEB−独立态反应自由能差分别为 `+2.0403335 eV` 与 `−1.5814870 eV`，不能解释为已验证的共吸附修正或用于移动 TS。详见 `NEB_MAPPING_VALIDATION.md`。

## 可复用的 Ti–Co 覆盖度根

已有 `zero_inventory_tico_diagnostic` 保存了一个受两个库存约束的高精度数值稳态候选：

- θ(1-己烯*) = `0.98953452226784795…`
- θ(OOH*) = `3.1966965438196039…×10⁻²⁰`
- θ(*) = `0.01046547773215204…`
- 条件性 C₆H₁₂O 净 TOF = `1.3236777105111254…×10⁻³⁶ s⁻¹`
- I₁=`2.2634×10⁻⁵⁷¹`，I₂=`0`；位点和 `1.0`。

高精度根复现与数值质量门控通过；但它仍是旧八步输入下的 provisional 根。局部稳定不证明从洁净表面全局可达；NEB 结构和热修正映射仍未通过物理核验。因此覆盖度筛选仅用于局部、探索性展示，不构成正式 TOF/OOH–TOF 结论。
