# `lyq_alkene_epoxidation` 代码逐行解读：Python 基础 + 项目作用

> 面向当前 `main` 分支。本文的阅读顺序是：**先理解这一行体现的 Python 基础概念，再理解它在 CatMAP/微观动力学项目中的实际作用**。这样既能复习 Python，又能快速建立对代码结构的认识。
>
> 说明：注释行也解释；**纯空行只起排版作用，不单独列出**。`.ipynb` 的 JSON 元数据、历史输出和 traceback 不属于实际执行源码，因此只解释代码单元中的 Python 语句。

## 1. 先看整体逻辑

```text
DFT/自由能原始数据
    ↓
data_only_not_executable.py / generate_input1.py
    ↓
形成能 E_form
    ↓
tini_energies.txt / tini_energies.csv
    ↓
多表面 energies.txt
    ↓
alkene_epoxidation.mkm
    ↓
CatMAP ReactionModel
    ↓
速率、覆盖度、生产速率、控制系数
    ↓
test.ipynb 后处理与绘图
```

### 建议的学习方法

看到每一行时，先问两个问题：

1. **Python 在做什么？**——导入、赋值、字典、循环、条件、函数、异常、文件 IO、对象属性……
2. **项目上为什么要这么做？**——形成能计算、CatMAP 命名映射、输入文件生成、模型配置、求解或后处理。

---

# 2. `1-generating_input_file/data_only_not_executable.py`

## 2.1 文件作用

这是一个“**单表面能量数据 → 形成能 → CatMAP TXT → 解析自检**”脚本。

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `from ase.symbols import string2symbols` | `from ... import ...` 表示从模块中只导入指定对象；之后可直接写 `string2symbols()`。 | 导入 ASE 的化学式解析函数，把 `H2O2`、`C6H12O` 等拆成原子序列，供形成能计算。 |
| 3 | `# abinitio_energies：` | `#` 后面的内容是**注释**，Python 不执行，用于给人阅读。 | 提示下面开始定义第一性原理能量数据。 |
| 4 | `# 1. 绝对 DFT 能量（不考虑温度，后续可加入频率修正）` | 注释可以记录变量的物理含义和假设。 | 说明输入可以是未做热修正的 DFT 总能。 |
| 5 | `# 2. 或已经得到的吉布斯自由能（后续无需再加频率）` | 注释不会改变变量值。 | 说明另一种输入可能已经是自由能，此时不应重复加热力学修正。 |
| 6 | `# 单位：eV` | 注释常用于记录单位，避免“数值正确但单位错误”。 | 指定所有能量按 eV 使用。 |
| 7 | `abinitio_energies = {` | `变量 = 值` 是**赋值**；`{}` 创建**字典 dict**，字典用“键:值”保存映射关系。 | 建立“物种 key → DFT/自由能”的总能数据库。 |
| 8 | `# ========== 气相物种 ==========` | 缩进内也可以写注释；不影响字典结构。 | 把下面几条数据标记为气相/储库物种。 |
| 9 | `'H2O2_gas': -18.122091,` | 字典元素格式是 `key: value`；字符串作 key，浮点数作 value；末尾逗号分隔元素。 | 保存 H2O2 气相能量。`gas` 后缀之后用于判断它不需要减裸板能。 |
| 10 | `'C6H12_gas': -95.031223,` | 同一字典可以连续保存多个同类型键值对。 | 保存 1-己烯气相能量。 |
| 11 | `'C6H12O_gas': -101.037771,` | 负号属于数值字面量的一部分。 | 保存环氧产物储库能量。 |
| 12 | `'H2O_gas': -14.204457,` | 字符串 key 必须唯一；重复 key 会覆盖前值。 | 保存水储库能量。 |
| 14 | `# ========== 吸附物种 ==========` | 注释用于人为分区。 | 下面开始保存“表面 + 吸附物”的总能。 |
| 15 | `'H2O2_111': -881.680784,` | 字典 key 本身可以编码结构化信息；这里用 `_` 人工约定“名称_site”。 | 保存 H2O2 在 111 表面的总能。 |
| 16 | `'Ha_111': -868.06986,` | Python 不理解 `Ha` 的化学意义，它只是普通字符串。 | 保存一个特殊 H 位点 `Ha` 的吸附体系总能。 |
| 17 | `'Hb_111': -867.91225,` | 同上，字符串名称的语义由程序员定义。 | 保存另一局域 H 位点 `Hb` 的总能。 |
| 18 | `'OOH_111': -877.1433925,` | 浮点数可保存较多小数位，但后续是否保留由 `round()` 决定。 | 保存 OOH* 中间体总能。 |
| 19 | `'C6H12_111': -958.3287775,` | 字典适合做“查表”，后续可通过 key 快速访问。 | 保存吸附己烯总能。 |
| 20 | `'C6H12O_111': -964.45279,` | key 与 value 没有限定类型，但本程序约定为 `str → float`。 | 保存吸附环氧产物总能。 |
| 21 | `'O_111': -866.997122,` | 字典中的元素顺序在现代 Python 中保持插入顺序，但逻辑不应依赖它。 | 保存 O* 总能。 |
| 22 | `'OH_111': -872.6800805,` | 这里仍是普通键值对。 | 保存 OH* 总能。 |
| 23 | `'H2O_111': -877.747855,` | 这里仍是普通键值对。 | 保存 H2O* 总能。 |
| 25 | `# ========== 过渡态 ==========` | 注释可以给数据增加语义分组。 | 下面两条是显式过渡态能量。 |
| 26 | `'OOH-C6H12_111': -983.369042,` | 字符串中可以包含 `-`；Python 不会自动把它当减号，因为它位于引号内。 | 保存环氧化过渡态总能。 |
| 27 | `'Hb-O_111': -872.064286,` | 同样，`Hb-O` 是字符串名称，不是 Python 运算式。 | 保存 H 转移过渡态总能。 |
| 29 | `# ========== 裸板 ==========` | 注释。 | 表明下面数据是干净表面的参考能。 |
| 30 | `'slab_111': -863.25438,` | 字典查值可以用 `dict[key]`。 | 裸板能后面既放入参考字典，又从吸附体系总能中扣除。 |
| 31 | `}` | `}` 结束字典字面量。 | 完成全部原始能量数据定义。 |
| 34 | `# 原子参考能，单位：eV` | 注释。 | 说明下面是形成能计算用的参考态。 |
| 35 | `ref_dict = {` | 再创建一个独立字典；变量名只要合法即可。 | 建立“元素/表面 → 参考能”的映射。 |
| 36 | `'O': -4.37,` | 字典 key 可用元素符号字符串。 | O 的原子参考能。 |
| 37 | `'C': -9.28,` | 同上。 | C 的原子参考能。 |
| 38 | `'H': -1.11,` | 同上。 | H 的原子参考能。 |
| 39 | `'111': abinitio_energies['slab_111'],` | `dict['key']` 是**按键索引字典**；右侧先求值，再存入新字典。 | 把裸板能直接作为 `111` 表面的参考能，避免手抄两遍。 |
| 40 | `}` | 结束字典。 | 完成参考能数据库。 |
| 42 | `# Ha/Hb 是不同吸附位置上的 H，不是标准化学元素符号。` | 注释。 | 解释为什么后面不能直接把 `Ha`/`Hb` 交给 ASE 解析。 |
| 43 | `# 因此显式指定这些特殊物种对应的真实化学组成。` | 注释。 | 引出特殊物种映射。 |
| 44 | `formula_map = {` | 字典也可做**名称转换表/映射表**。 | 把项目中的非标准名字映射到真正的化学式。 |
| 45 | `'Ha': 'H',` | 字典 value 也可以是字符串。 | `Ha` 实际就是一个 H。 |
| 46 | `'Hb': 'H',` | 同上。 | `Hb` 实际也是一个 H。 |
| 47 | `'Hb-O': 'HO',` | key/value 都是普通字符串。 | `Hb-O` 在元素组成上是 H + O。 |
| 48 | `}` | 结束字典。 | 完成特殊名称映射。 |
| 51 | `def get_formation_energies(energy_dict, references):` | `def` **定义函数**；括号中是形参。函数只有被调用时才真正执行函数体。 | 定义通用形成能计算器。 |
| 52 | `"""根据裸板能和原子参考能计算形成能。"""` | 函数体第一条三引号字符串是**docstring 文档字符串**，可通过 `.__doc__` 查看。 | 给函数说明用途。 |
| 53 | `formation_energies = {}` | `{}` 创建空字典；变量位于函数内部，是局部变量。 | 准备保存计算结果。 |
| 55 | `for key, energy in energy_dict.items():` | `for` 是循环；`dict.items()` 返回键值对；`key, energy` 是**序列解包**。 | 逐个处理所有物种及其原始能量。 |
| 56 | `if '_' not in key:` | `if` 做条件判断；`not in` 检查元素/子串是否不存在。 | 验证 key 是否符合“物种_site”命名规则。 |
| 57 | `raise ValueError(f'物种 key 缺少 site 信息: {key}')` | `raise` 主动抛异常；`ValueError` 表示值格式不符合预期；`f'...{变量}...'` 是 **f-string**。 | 一旦命名非法就立即停止，并把问题 key 打进错误信息。 |
| 59 | `# 只按第一个下划线拆分，避免物种名中出现其他字符时拆分过多。` | 注释。 | 解释此前 `too many values to unpack` 类问题的修复思路。 |
| 60 | `name, site = key.split('_', 1)` | `str.split(sep, maxsplit)` 拆字符串；`maxsplit=1` 最多切一次；左边两个变量做解包。 | 安全得到物种名和 site，避免 key 中更多下划线时拆成 3 段以上。 |
| 62 | `if 'slab' in name:` | `in` 对字符串可检查子串。 | 判断当前条目是不是裸板参考。 |
| 63 | `continue` | `continue` 立即结束本轮循环，直接进入下一轮。 | 裸板只做参考，不生成形成能物种条目。 |
| 65 | `E0 = energy` | 简单赋值；数字是不可变对象，这里可把 `E0` 当工作变量。 | 从原始能量开始逐步扣除参考能。 |
| 67 | `# 表面吸附态和过渡态减去裸板能，气相不减。` | 注释。 | 说明形成能公式的分支逻辑。 |
| 68 | `if site == '111':` | `==` 判断两值是否相等。 | 只对表面态/过渡态执行裸板扣除。 |
| 69 | `E0 -= references['111']` | `-=` 是**增强赋值**，等价于 `E0 = E0 - ...`。 | 从“表面+物种”总能中减去干净表面能。 |
| 71 | `formula = formula_map.get(name, name.replace('-', ''))` | `dict.get(key, default)` 在 key 不存在时返回默认值；`str.replace()` 返回替换后的新字符串。 | 特殊名用映射表；普通过渡态名去掉 `-` 后作为化学式。 |
| 73 | `try:` | `try/except` 是**异常处理**结构，用来捕获可能失败的代码。 | 准备捕获 ASE 化学式解析失败。 |
| 74 | `composition = string2symbols(formula)` | 函数调用：把实参传给函数并接收返回值。 | 例如把 `H2O2` 解析成 H、H、O、O 序列。 |
| 75 | `except Exception as exc:` | `except` 捕获异常；`as exc` 把异常对象保存到变量。`Exception` 是较宽泛的异常基类。 | 任何化学式解析问题都进入统一错误提示。 |
| 76 | `raise ValueError(` | 函数调用可以跨多行书写，只要括号没有闭合。 | 开始重新包装异常。 |
| 77 | `f'无法解析物种 {name}，对应化学式 {formula}，原 key: {key}'` | f-string 可在字符串中直接嵌变量表达式。 | 错误信息同时给出项目名、化学式和原 key。 |
| 78 | `) from exc` | `raise ... from exc` 建立**异常链**，保留原始异常作为根因。 | 既有友好提示，又保留 ASE 的底层 traceback。 |
| 80 | `for atom in composition:` | `for` 可直接遍历列表中的每个元素。 | 对化学式中的每个原子逐个扣参考能；同元素出现多次就重复扣。 |
| 81 | `if atom not in references:` | `not in` 对字典默认检查 key。 | 确保每个原子都有参考能。 |
| 82 | `raise KeyError(f'物种 {key} 中的原子 {atom} 未在 ref_dict 中定义参考能')` | `KeyError` 通常用于“字典缺少需要的键”。 | 缺失参考能时直接阻止错误形成能继续传播。 |
| 83 | `E0 -= references[atom]` | 字典索引 + 增强赋值。 | 减去当前原子的参考能。 |
| 85 | `formation_energies[key] = round(E0, 3)` | `dict[key] = value` 写入/更新字典；`round(x, 3)` 保留 3 位小数。 | 把最终形成能存入结果，控制输出精度。 |
| 87 | `return formation_energies` | `return` 结束函数并把结果返回调用者。 | 返回完整形成能字典。 |
| 90 | `formation_energies = get_formation_energies(abinitio_energies, ref_dict)` | 调用函数时，实参与形参按位置对应；返回值再赋给变量。 | 真正执行形成能计算。 |
| 92 | `for key, value in formation_energies.items():` | 字典遍历 + 两变量解包。 | 遍历所有计算结果用于检查。 |
| 93 | `print(f'{key} {value}')` | `print()` 输出到终端；f-string 组织文本。 | 打印物种及形成能，便于人工核对。 |
| 96 | `frequency_dict = {` | 再创建一个字典。 | 建立“物种 → 振动频率列表”数据库。 |
| 97 | `'H2O2_gas': [],` | `[]` 是空列表 list。 | H2O2 暂无频率数据。 |
| 98 | `'C6H12_gas': [],` | 空列表是可变容器，可以后追加元素。 | C6H12 暂无频率。 |
| 99 | `'C6H12O_gas': [],` | 同上。 | C6H12O 暂无频率。 |
| 100 | `'H2O_gas': [],` | 同上。 | H2O 暂无频率。 |
| 101 | `'H2O2_111': [],` | 同上。 | H2O2* 暂无频率。 |
| 102 | `'Ha_111': [],` | 同上。 | Ha* 暂无频率。 |
| 103 | `'Hb_111': [],` | 同上。 | Hb* 暂无频率。 |
| 104 | `'OOH_111': [],` | 同上。 | OOH* 暂无频率。 |
| 105 | `'C6H12_111': [],` | 同上。 | C6H12* 暂无频率。 |
| 106 | `'C6H12O_111': [],` | 同上。 | C6H12O* 暂无频率。 |
| 107 | `'O_111': [],` | 同上。 | O* 暂无频率。 |
| 108 | `'OH_111': [],` | 同上。 | OH* 暂无频率。 |
| 109 | `'H2O_111': [],` | 同上。 | H2O* 暂无频率。 |
| 110 | `'OOH-C6H12_111': [],` | 同上。 | 第一过渡态暂无频率。 |
| 111 | `'Hb-O_111': [],` | 同上。 | 第二过渡态暂无频率。 |
| 112 | `'slab_111': [],` | 同上。 | 裸板频率也暂留空。 |
| 113 | `}` | 结束字典。 | 完成频率数据库。 |
| 116 | `def make_input_file(file_name, energy_dict, frequencies):` | 定义带 3 个参数的函数。 | 封装 CatMAP TXT 输入文件生成逻辑。 |
| 117 | `"""生成 CatMAP TableParser 使用的制表符分隔输入文件。"""` | docstring。 | 说明函数输出格式。 |
| 118 | `header = '\t'.join([` | `str.join(iterable)` 用指定分隔符连接多个字符串；`\t` 是制表符。 | 开始生成 CatMAP 表头。 |
| 119 | `'surface_name',` | 列表中的字符串元素。 | 第一列：表面名称。 |
| 120 | `'site_name',` | 同上。 | 第二列：位点名称。 |
| 121 | `'species_name',` | 同上。 | 第三列：物种名。 |
| 122 | `'formation_energy',` | 同上。 | 第四列：形成能。 |
| 123 | `'frequencies',` | 同上。 | 第五列：频率。 |
| 124 | `'reference',` | 同上。 | 第六列：数据来源。 |
| 125 | `])` | `]` 结束列表，`)` 结束 `join()` 调用。 | 得到最终表头字符串。 |
| 127 | `lines = []` | 创建空列表。 | 准备收集每一个物种的数据行。 |
| 129 | `for key, energy in energy_dict.items():` | 字典遍历 + 解包。 | 逐物种生成输出行。 |
| 130 | `if '_' not in key:` | 条件判断 + 成员测试。 | 非法 key 不参与输出。 |
| 131 | `continue` | 跳过当前循环。 | 直接处理下一个物种。 |
| 133 | `# 必须只拆分一次；不能使用 key.split('_')。` | 注释。 | 强调命名拆分规则。 |
| 134 | `name, site = key.split('_', 1)` | 字符串拆分 + 解包。 | 得到 species_name 和 site_name。 |
| 136 | `if 'slab' in name:` | 条件判断。 | 识别裸板。 |
| 137 | `continue` | 跳过当前循环。 | 裸板不写入 CatMAP 物种表。 |
| 139 | `frequency = frequencies.get(key, [])` | `dict.get()` 可避免 key 不存在时报 `KeyError`。 | 若该物种没频率，就用空列表。 |
| 140 | `surface = None if site == 'gas' else 'tini'` | 这是**条件表达式/三元表达式**：`A if 条件 else B`。`None` 表示“无值/空对象”。 | 气相没有表面；非气相统一写 `tini`。这里也是当前脚本的硬编码点。 |
| 141 | `outline = [` | 创建列表。 | 开始组织当前物种的一整行字段。 |
| 142 | `surface,` | 变量作为列表元素。 | surface_name。 |
| 143 | `site,` | 同上。 | site_name。 |
| 144 | `name,` | 同上。 | species_name。 |
| 145 | `energy,` | 同上。 | formation_energy。 |
| 146 | `frequency,` | 同上。 | frequencies。 |
| 147 | `'Input File Tutorial.',` | 字符串常量。 | 固定 reference 说明。 |
| 148 | `]` | 结束列表。 | 当前行字段组织完成。 |
| 149 | `lines.append('\t'.join(str(value) for value in outline))` | `append()` 往列表末尾加一个元素；`str(value) for value in outline` 是**生成器表达式**；`join()` 再连接。 | 把所有字段转字符串并按制表符拼成一行。 |
| 151 | `lines.sort()` | `list.sort()` 就地排序列表。 | 让输出顺序稳定，方便比较不同运行结果。 |
| 152 | `input_file = '\n'.join([header] + lines)` | `+` 可拼接两个列表；`\n` 是换行符。 | 把表头和所有数据行拼成完整文本。 |
| 154 | `with open(file_name, 'w', encoding='utf-8') as file:` | `open()` 打开文件；`'w'` 是覆盖写入；`with` 是**上下文管理器**，退出代码块会自动关闭文件。 | 安全创建 CatMAP TXT 文件。 |
| 155 | `file.write(input_file)` | 文件对象的 `write()` 写字符串。 | 把完整能量表写入磁盘。 |
| 157 | `print(f'Successfully created input file: {file_name}')` | `print()` + f-string。 | 告知用户文件生成成功。 |
| 160 | `file_name = 'tini_energies.txt'` | 变量赋值。 | 指定实际输出文件名。 |
| 161 | `make_input_file(file_name, formation_energies, frequency_dict)` | 函数调用。 | 真正生成 `tini_energies.txt`。 |
| 164 | `# 测试 CatMAP 是否能够正确解析生成的输入文件。` | 注释。 | 下面进入“输出文件自检”。 |
| 165 | `from catmap.model import ReactionModel` | 从子模块导入类。 | 获取 CatMAP 核心模型类。 |
| 166 | `from catmap.parsers import TableParser` | 从模块导入类。 | 获取表格输入解析器。 |
| 168 | `rxm = ReactionModel()` | `类名()` 实例化对象；变量保存对象引用。 | 创建最小 CatMAP 模型用于测试。 |
| 169 | `rxm.surface_names = ['tini']` | `对象.属性 = 值` 是**属性赋值**。 | 声明模型表面为 tini。 |
| 171 | `# CatMAP 内部表面物种 key 使用 *_s 形式；输入表中的 species_name 仍保持` | 注释。 | 解释 CatMAP 内部命名与输入表命名的差别。 |
| 172 | `# H2O2、Ha、Hb 等名称，由 TableParser 根据 site_name=111 匹配到 s 位点。` | 注释。 | 说明 `_s` 与 `111` 之间的映射关系。 |
| 173 | `rxm.adsorbate_names = (` | `()` 创建元组；多行书写可提升可读性。 | 开始声明所有吸附物种。 |
| 174 | `'H2O2_s',` | 元组元素；单个字符串后面的逗号很重要。 | 声明 H2O2*。 |
| 175 | `'Ha_s',` | 同上。 | 声明 Ha*。 |
| 176 | `'Hb_s',` | 同上。 | 声明 Hb*。 |
| 177 | `'OOH_s',` | 同上。 | 声明 OOH*。 |
| 178 | `'C6H12_s',` | 同上。 | 声明 C6H12*。 |
| 179 | `'C6H12O_s',` | 同上。 | 声明 C6H12O*。 |
| 180 | `'O_s',` | 同上。 | 声明 O*。 |
| 181 | `'OH_s',` | 同上。 | 声明 OH*。 |
| 182 | `'H2O_s',` | 同上。 | 声明 H2O*。 |
| 183 | `)` | 结束元组。 | adsorbate_names 定义完成。 |
| 184 | `rxm.transition_state_names = ('OOH-C6H12_s', 'Hb-O_s')` | 元组 + 对象属性赋值。 | 声明两个显式过渡态。 |
| 185 | `rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')` | 元组 + 属性赋值。 | 声明储库/气相物种。 |
| 186 | `rxm.site_names = ('s',)` | **单元素元组必须写逗号**，`('s')` 只是字符串，`('s',)` 才是 tuple。 | 声明模型只有一个抽象位点 `s`。 |
| 188 | `# 非标准名称必须显式定义元素组成。` | 注释。 | 引出特殊物种配置。 |
| 189 | `rxm.species_definitions = {` | 给对象属性赋一个嵌套字典。 | 手工补充位点映射和元素组成。 |
| 190 | `'s': {'site_names': ['111']},` | **嵌套数据结构**：字典中放字典，内层还放列表。 | 把抽象 `s` 位点映射到输入表中的 `111`。 |
| 191 | `'Ha_s': {'composition': {'H': 1}},` | 多层嵌套字典。 | 指定 Ha_s = H1。 |
| 192 | `'Hb_s': {'composition': {'H': 1}},` | 同上。 | 指定 Hb_s = H1。 |
| 193 | `'Hb-O_s': {'composition': {'H': 1, 'O': 1}},` | 同一内层字典可有多个键。 | 指定 Hb-O_s = H1O1。 |
| 194 | `}` | 结束外层字典。 | 完成 species_definitions。 |
| 196 | `parser = TableParser(rxm)` | 实例化类并把已有对象作为参数传入。 | 创建与 `rxm` 绑定的表格解析器。 |
| 197 | `parser.input_file = file_name` | 对象属性赋值。 | 告诉解析器读取刚生成的 TXT。 |
| 198 | `parser.parse()` | 调用对象方法。 | 真正解析输入文件；列名/物种/能量不匹配会在此报错。 |
| 200 | `print('CatMAP parsing successful.')` | `print()` 输出固定字符串。 | 只有解析成功才会执行到这里。 |
| 201 | `for key, value in rxm.species_definitions.items():` | 遍历对象内部字典。 | 检查解析后 CatMAP 如何保存各物种信息。 |
| 202 | `print(f'{key} {value}')` | f-string + 输出。 | 把解析结果打印出来。 |

---

# 3. `1-generating_input_file/generate_input1.py`

## 3.1 文件作用

这是上一个脚本的扩展版：增加 **CSV 输出**，并用 `try/except` 对 CatMAP 解析失败进行更清楚的报错。

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `# 生成 CatMAP 输入文件：计算形成能、输出 TXT/CSV，并验证解析结果。` | 注释。 | 用一句话概括脚本职责。 |
| 2 | `from ase.symbols import string2symbols` | 指定对象导入。 | 用 ASE 解析化学式。 |
| 3 | `import csv` | `import 模块` 导入整个模块；后续通过 `csv.xxx` 使用。`csv` 属于 Python 标准库。 | 用标准库生成规范 CSV。 |
| 5 | `abinitio_energies = {` | 字典定义。 | 开始保存该表面的原始能量。 |
| 6 | `# ========== 气相物种 ==========` | 注释。 | 标记气相物种。 |
| 7 | `'H2O2_gas': -18.122091,` | 字典键值对。 | H2O2 气相能。 |
| 8 | `'C6H12_gas': -95.031223,` | 字典键值对。 | C6H12 气相能。 |
| 9 | `'C6H12O_gas': -101.037771,` | 字典键值对。 | C6H12O 气相能。 |
| 10 | `'H2O_gas': -14.204457,` | 字典键值对。 | H2O 气相能。 |
| 12 | `# ========== 吸附物种 ==========` | 注释。 | 标记吸附体系。 |
| 13 | `'H2O2_111': -891.610975,` | 字典键值对。 | H2O2* 总能。 |
| 14 | `'Ha_111': -874.30923,` | 字典键值对。 | Ha* 总能。 |
| 15 | `'Hb_111': -879.2156865,` | 字典键值对。 | Hb* 总能。 |
| 16 | `'OOH_111': -887.017207,` | 字典键值对。 | OOH* 总能。 |
| 17 | `'C6H12_111': -968.7609275,` | 字典键值对。 | C6H12* 总能。 |
| 18 | `'C6H12O_111': -974.959544,` | 字典键值对。 | C6H12O* 总能。 |
| 19 | `'O_111': -876.0788065,` | 字典键值对。 | O* 总能。 |
| 20 | `'OH_111': -881.650587,` | 字典键值对。 | OH* 总能。 |
| 21 | `'H2O_111': -887.646684,` | 字典键值对。 | H2O* 总能。 |
| 23 | `# ========== 过渡态 ==========` | 注释。 | 标记过渡态。 |
| 24 | `'OOH-C6H12_111': -983.35808,` | 引号中的 `-` 只是字符串字符。 | 环氧化 TS 总能。 |
| 25 | `'Hb-O_111': -881.739798,` | 字典键值对。 | H 转移 TS 总能。 |
| 27 | `# ========== 裸板 ==========` | 注释。 | 标记裸板。 |
| 28 | `'slab_111': -873.630641,` | 字典键值对。 | 保存裸板参考能。 |
| 29 | `}` | 结束字典。 | 原始能量输入结束。 |
| 31 | `# 预计算原子参考能，单位：eV` | 注释。 | 说明参考能。 |
| 32 | `ref_dict = {` | 字典定义。 | 建立形成能参考态。 |
| 33 | `'O': -4.37,` | 字典键值对。 | O 参考能。 |
| 34 | `'C': -9.28,` | 字典键值对。 | C 参考能。 |
| 35 | `'H': -1.11,` | 字典键值对。 | H 参考能。 |
| 36 | `'111': abinitio_energies['slab_111'],` | 字典按 key 读取值。 | 直接复用裸板能作为表面参考。 |
| 37 | `}` | 结束字典。 | 参考字典完成。 |
| 39 | `# 特殊物种名称 -> 实际化学组成。` | 注释。 | 引出名称映射。 |
| 40 | `# Ha/Hb 表示不同吸附位置上的 H，不是化学元素 Ha/Hb。` | 注释。 | 说明非标准命名问题。 |
| 41 | `formula_map = {` | 字典作为映射表。 | 保存特殊物种真实化学式。 |
| 42 | `'Ha': 'H',` | 字符串到字符串映射。 | Ha → H。 |
| 43 | `'Hb': 'H',` | 同上。 | Hb → H。 |
| 44 | `'Hb-O': 'HO',` | 同上。 | Hb-O → HO。 |
| 45 | `}` | 结束字典。 | 映射完成。 |
| 48 | `def get_formation_energies(energy_dict, references):` | 函数定义。 | 定义形成能计算器。 |
| 49 | `"""根据原子参考能和裸板能计算形成能。"""` | docstring。 | 函数说明。 |
| 50 | `formation_energies = {}` | 创建局部空字典。 | 保存结果。 |
| 52 | `for key, energy in energy_dict.items():` | 字典遍历 + 解包。 | 遍历全部物种。 |
| 53 | `if '_' not in key:` | `if` + `not in`。 | 检查命名格式。 |
| 54 | `raise ValueError(f'物种 key 缺少 site 信息: {key}')` | 主动抛异常 + f-string。 | 非法 key 立即报错。 |
| 56 | `name, site = key.split('_', 1)` | `split()` + 多变量解包。 | 只切一次，避免多下划线拆分错误。 |
| 58 | `if 'slab' in name:` | 成员判断。 | 识别裸板。 |
| 59 | `continue` | 循环控制。 | 裸板不进入形成能结果。 |
| 61 | `E0 = energy` | 变量赋值。 | 建立可逐步修改的工作能量。 |
| 63 | `# 表面吸附态/过渡态先减去裸板能；气相不减裸板。` | 注释。 | 说明物理公式。 |
| 64 | `if site == '111':` | 相等比较。 | 仅表面物种减裸板。 |
| 65 | `E0 -= references['111']` | 增强赋值。 | 扣除裸板能。 |
| 67 | `# 特殊名称显式映射，普通名称则移除 TS 连字符后按化学式解析。` | 注释。 | 说明名称处理策略。 |
| 68 | `formula = formula_map.get(name, name.replace('-', ''))` | `dict.get()` + 字符串替换。 | 得到 ASE 可解析的化学式。 |
| 70 | `try:` | 异常捕获结构开始。 | 捕获解析失败。 |
| 71 | `composition = string2symbols(formula)` | 函数调用。 | 解析元素组成。 |
| 72 | `except Exception as exc:` | 捕获异常对象。 | 接管 ASE 解析错误。 |
| 73 | `raise ValueError(` | 跨行构造异常。 | 改成更直观的项目级错误。 |
| 74 | `f'无法解析物种 {name}，对应化学式 {formula}，原 key: {key}'` | f-string。 | 把关键变量写入错误信息。 |
| 75 | `) from exc` | 异常链。 | 保留原始原因。 |
| 77 | `for atom in composition:` | 遍历序列。 | 逐原子处理。 |
| 78 | `if atom not in references:` | 字典 key 成员检查。 | 验证参考能完整。 |
| 79 | `raise KeyError(f'物种 {key} 中的原子 {atom} 未在 ref_dict 中定义参考能')` | `KeyError` + f-string。 | 缺参考能时停止。 |
| 80 | `E0 -= references[atom]` | 字典索引 + `-=`。 | 扣元素参考能。 |
| 82 | `formation_energies[key] = round(E0, 3)` | 字典写值 + `round()`。 | 保存 3 位小数形成能。 |
| 84 | `return formation_energies` | 返回函数结果。 | 把结果传出。 |
| 87 | `formation_energies = get_formation_energies(abinitio_energies, ref_dict)` | 函数调用并接收返回值。 | 执行形成能计算。 |
| 89 | `# 打印形成能检查` | 注释。 | 标记人工检查阶段。 |
| 90 | `for key, value in formation_energies.items():` | 遍历字典。 | 查看所有结果。 |
| 91 | `print(f'{key} {value}')` | 输出 + f-string。 | 打印形成能。 |
| 93 | `frequency_dict = {` | 字典。 | 建立频率表。 |
| 94 | `'H2O2_gas': [],` | 空列表。 | 暂无频率。 |
| 95 | `'C6H12_gas': [],` | 空列表。 | 暂无频率。 |
| 96 | `'C6H12O_gas': [],` | 空列表。 | 暂无频率。 |
| 97 | `'H2O_gas': [],` | 空列表。 | 暂无频率。 |
| 98 | `'H2O2_111': [],` | 空列表。 | 暂无频率。 |
| 99 | `'Ha_111': [],` | 空列表。 | 暂无频率。 |
| 100 | `'Hb_111': [],` | 空列表。 | 暂无频率。 |
| 101 | `'OOH_111': [],` | 空列表。 | 暂无频率。 |
| 102 | `'C6H12_111': [],` | 空列表。 | 暂无频率。 |
| 103 | `'C6H12O_111': [],` | 空列表。 | 暂无频率。 |
| 104 | `'O_111': [],` | 空列表。 | 暂无频率。 |
| 105 | `'OH_111': [],` | 空列表。 | 暂无频率。 |
| 106 | `'H2O_111': [],` | 空列表。 | 暂无频率。 |
| 107 | `'OOH-C6H12_111': [],` | 空列表。 | 暂无频率。 |
| 108 | `'Hb-O_111': [],` | 空列表。 | 暂无频率。 |
| 109 | `'slab_111': [],` | 空列表。 | 暂无频率。 |
| 110 | `}` | 结束字典。 | 频率表完成。 |
| 113 | `def make_input_file(file_name, energy_dict, frequencies):` | 函数定义。 | 定义 TXT 生成器。 |
| 114 | `"""生成 CatMAP TableParser 使用的制表符分隔输入文件。"""` | docstring。 | 说明函数功能。 |
| 115 | `header = '\t'.join([` | `join()` + 列表。 | 开始构造表头。 |
| 116 | `'surface_name',` | 列表元素。 | 表面列。 |
| 117 | `'site_name',` | 列表元素。 | 位点列。 |
| 118 | `'species_name',` | 列表元素。 | 物种列。 |
| 119 | `'formation_energy',` | 列表元素。 | 形成能列。 |
| 120 | `'frequencies',` | 列表元素。 | 频率列。 |
| 121 | `'reference',` | 列表元素。 | 来源列。 |
| 122 | `])` | 结束列表和函数调用。 | 表头完成。 |
| 124 | `lines = []` | 空列表。 | 保存输出行。 |
| 126 | `for key, energy in energy_dict.items():` | 字典循环。 | 逐物种输出。 |
| 127 | `if '_' not in key:` | 条件判断。 | 过滤非法 key。 |
| 128 | `continue` | 跳过本轮。 | 不输出非法条目。 |
| 130 | `name, site = key.split('_', 1)` | 字符串拆分。 | 得到名称和 site。 |
| 132 | `if 'slab' in name:` | 条件判断。 | 识别裸板。 |
| 133 | `continue` | 跳过本轮。 | 不输出裸板。 |
| 135 | `frequency = frequencies.get(key, [])` | `dict.get()` 默认值。 | 频率不存在就用空列表。 |
| 136 | `surface = None if site == 'gas' else 'tini'` | 三元表达式。 | 气相 surface=None；表面固定 tini。 |
| 137 | `outline = [surface, site, name, energy, frequency, 'Input File Tutorial.']` | 一行即可创建列表，元素可以是不同类型。 | 把 6 个输出字段组织成一行。 |
| 138 | `lines.append('\t'.join(str(value) for value in outline))` | `append()` + 生成器表达式 + `join()`。 | 转为 CatMAP 制表符格式并加入结果。 |
| 140 | `lines.sort()` | 列表就地排序。 | 输出顺序稳定。 |
| 141 | `input_file = '\n'.join([header] + lines)` | 列表拼接 + 字符串连接。 | 生成完整 TXT 内容。 |
| 143 | `with open(file_name, 'w', encoding='utf-8') as file:` | 上下文管理器 + 文件写入模式。 | 创建 TXT。 |
| 144 | `file.write(input_file)` | 文件方法调用。 | 写入数据。 |
| 146 | `print(f'Successfully created input file: {file_name}')` | 输出。 | 打印成功提示。 |
| 149 | `def make_csv_file(file_name, energy_dict, frequencies):` | 函数定义。 | 定义 CSV 生成器。 |
| 150 | `"""将同一套形成能数据输出为 CSV，便于表格软件检查。"""` | docstring。 | 说明 CSV 用途。 |
| 151 | `fieldnames = [` | 列表定义。 | 准备 CSV 列名。 |
| 152 | `'surface_name',` | 列表元素。 | CSV 表面列。 |
| 153 | `'site_name',` | 列表元素。 | CSV 位点列。 |
| 154 | `'species_name',` | 列表元素。 | CSV 物种列。 |
| 155 | `'formation_energy',` | 列表元素。 | CSV 能量列。 |
| 156 | `'frequencies',` | 列表元素。 | CSV 频率列。 |
| 157 | `'reference',` | 列表元素。 | CSV 来源列。 |
| 158 | `]` | 结束列表。 | 列名完成。 |
| 160 | `with open(file_name, 'w', newline='', encoding='utf-8') as csvfile:` | 文件上下文管理器；写 CSV 时常用 `newline=''` 避免 Windows 多空行。 | 创建标准 CSV 文件。 |
| 161 | `writer = csv.DictWriter(csvfile, fieldnames=fieldnames)` | `模块.类()` 实例化对象；关键字参数 `fieldnames=...` 提高可读性。 | 创建按“字段名字典”写 CSV 的 writer。 |
| 162 | `writer.writeheader()` | 调用对象方法。 | 写入 CSV 表头。 |
| 164 | `for key, energy in energy_dict.items():` | 字典遍历。 | 逐物种写 CSV。 |
| 165 | `if '_' not in key:` | 条件判断。 | 过滤非法 key。 |
| 166 | `continue` | 跳过本轮。 | 不写非法条目。 |
| 168 | `name, site = key.split('_', 1)` | 字符串拆分。 | 提取名称和 site。 |
| 170 | `if 'slab' in name:` | 条件判断。 | 识别裸板。 |
| 171 | `continue` | 跳过本轮。 | 裸板不写入物种表。 |
| 173 | `freq_list = frequencies.get(key, [])` | `dict.get()`。 | 获取频率列表。 |
| 174 | `freq_str = ';'.join(` | `join()` 可把一系列字符串合成单字符串。 | 准备把多个频率塞进一个 CSV 单元格。 |
| 175 | `str(int(freq))` | 类型转换可嵌套：先 `int()` 再 `str()`。 | 对整数频率去掉多余 `.0`。 |
| 176 | `if isinstance(freq, (int, float)) and float(freq).is_integer()` | `isinstance()` 做类型判断；`and` 是逻辑与；`.is_integer()` 判断浮点值是否为整数。 | 决定一个频率是否可安全按整数输出。 |
| 177 | `else str(freq)` | 条件表达式的 `else` 分支。 | 非整数频率保留原表示。 |
| 178 | `for freq in freq_list` | 这是生成器表达式中的 `for`。 | 对每个频率执行同一格式化规则。 |
| 179 | `) if freq_list else ''` | 外层又使用条件表达式；空列表在布尔上下文中为 False。 | 没有频率时 CSV 单元格写空字符串。 |
| 181 | `writer.writerow({` | 传一个字典给方法调用。 | 开始写当前物种的一行 CSV。 |
| 182 | `'surface_name': None if site == 'gas' else 'tini',` | 字典 value 可以直接写条件表达式。 | 气相无表面，表面态写 tini。 |
| 183 | `'site_name': site,` | 字典键值对。 | 写 site。 |
| 184 | `'species_name': name,` | 字典键值对。 | 写物种名。 |
| 185 | `'formation_energy': energy,` | 字典键值对。 | 写形成能。 |
| 186 | `'frequencies': freq_str,` | 字典键值对。 | 写格式化频率。 |
| 187 | `'reference': 'Input File Tutorial.',` | 字典键值对。 | 写来源说明。 |
| 188 | `})` | 结束字典并结束 `writerow()` 调用。 | 当前 CSV 行写完。 |
| 190 | `print(f'Successfully created CSV file: {file_name}')` | 输出 + f-string。 | 提示 CSV 生成成功。 |
| 193 | `# 生成输出文件` | 注释。 | 下面真正调用生成函数。 |
| 194 | `txt_file_name = 'tini_energies.txt'` | 变量赋值。 | TXT 文件名。 |
| 195 | `csv_file_name = 'tini_energies.csv'` | 变量赋值。 | CSV 文件名。 |
| 196 | `make_input_file(txt_file_name, formation_energies, frequency_dict)` | 函数调用。 | 生成 TXT。 |
| 197 | `make_csv_file(csv_file_name, formation_energies, frequency_dict)` | 函数调用。 | 生成 CSV。 |
| 199 | `# CatMAP 解析测试` | 注释。 | 下面执行输入文件自检。 |
| 200 | `try:` | 异常处理开始。 | 把可能失败的 CatMAP 解析包起来。 |
| 201 | `from catmap.model import ReactionModel` | 导入可以放在代码块内部，只有执行到这里时才发生。 | 导入模型类。 |
| 202 | `from catmap.parsers import TableParser` | 同上。 | 导入解析器类。 |
| 204 | `rxm = ReactionModel()` | 类实例化。 | 创建空模型。 |
| 205 | `rxm.surface_names = ['tini']` | 对象属性赋值。 | 声明表面。 |
| 207 | `# CatMAP 内部物种 key 必须带 site 后缀。` | 注释。 | 解释 `_s`/`_g` 命名。 |
| 208 | `# 输入表中的 species_name 仍保持 H2O2/Ha/Hb/...，TableParser 会根据 site_name=111` | 注释。 | 说明表格中不直接写 `_s`。 |
| 209 | `# 将这些名称匹配到 *_s 物种。` | 注释。 | 说明 TableParser 的映射职责。 |
| 210 | `rxm.adsorbate_names = (` | 元组赋给对象属性。 | 开始声明吸附物集合。 |
| 211 | `'H2O2_s', 'Ha_s', 'Hb_s', 'OOH_s',` | 同一行可写多个元组元素。 | 前四个吸附物。 |
| 212 | `'C6H12_s', 'C6H12O_s', 'O_s', 'OH_s', 'H2O_s'` | 元组元素。 | 其余吸附物。 |
| 213 | `)` | 结束元组。 | 吸附物列表完成。 |
| 214 | `rxm.transition_state_names = ('OOH-C6H12_s', 'Hb-O_s')` | 元组 + 属性。 | 声明两个 TS。 |
| 215 | `rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')` | 元组 + 属性。 | 声明储库物种。 |
| 216 | `rxm.site_names = ('s',)` | 单元素元组语法。 | 声明一个 site 类型。 |
| 218 | `# Ha/Hb/Hb-O 不是标准化学式，必须在与模型一致的 *_s key 上显式指定组成。` | 注释。 | 解释特殊组成定义。 |
| 219 | `rxm.species_definitions = {` | 嵌套字典赋给对象属性。 | 开始补充物种元数据。 |
| 220 | `'s': {'site_names': ['111']},` | 嵌套字典/列表。 | s → 111 位点映射。 |
| 221 | `'Ha_s': {'composition': {'H': 1}},` | 多层嵌套字典。 | Ha_s = H1。 |
| 222 | `'Hb_s': {'composition': {'H': 1}},` | 多层嵌套字典。 | Hb_s = H1。 |
| 223 | `'Hb-O_s': {'composition': {'H': 1, 'O': 1}},` | 多层嵌套字典。 | Hb-O_s = H1O1。 |
| 224 | `}` | 结束字典。 | 配置完成。 |
| 226 | `parser = TableParser(rxm)` | 类实例化。 | 创建解析器。 |
| 227 | `parser.input_file = txt_file_name` | 属性赋值。 | 指定 TXT。 |
| 228 | `parser.parse()` | 方法调用。 | 真正解析。 |
| 230 | `print('CatMAP parsing successful.')` | 输出。 | 成功提示。 |
| 231 | `for key, value in rxm.species_definitions.items():` | 字典遍历。 | 查看解析结果。 |
| 232 | `print(f'{key} {value}')` | 输出。 | 打印物种定义。 |
| 234 | `except Exception as exc:` | 捕获 `try` 中异常并绑定给 `exc`。 | 统一处理 CatMAP 解析失败。 |
| 235 | `print('CatMAP 解析测试失败：')` | 输出。 | 给出中文失败标题。 |
| 236 | `print(exc)` | 异常对象可以直接 `print()`。 | 快速显示错误消息。 |
| 237 | `# 不再吞掉异常：解析失败时让进程返回非 0 退出码，便于定位问题。` | 注释。 | 说明为什么下一行还要再次抛异常。 |
| 238 | `raise` | 在 `except` 中单独写 `raise` 会**重新抛出当前异常**，保留原 traceback。 | 让脚本真正以失败结束，便于调试和自动化检测。 |

### 这一文件最值得掌握的 Python 基础

`dict`、`list`、函数、`for`、`if`、`split()`、`try/except/raise`、f-string、生成器表达式、`with open()`、`csv.DictWriter`、对象属性。

---

# 4. `2-creating_microkinetic_model/alkene_epoxidation.mkm`

`.mkm` 虽然扩展名不是 `.py`，但 CatMAP 把它当作 **Python 风格配置脚本**读取，所以里面的大部分语法就是 Python 的变量、列表、字典和字符串。

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `#` | 只有 `#` 的行也是注释。 | 视觉分隔。 |
| 2 | `# Microkinetic model parameters` | 注释。 | 标记模型参数区。 |
| 3 | `#` | 注释。 | 分隔。 |
| 5 | `# 1-hexene epoxidation with H2O2` | 注释。 | 指明反应体系。 |
| 6 | `rxn_expressions = [` | 列表赋值。 | 开始定义基元反应网络。 |
| 7 | `'* + H2O2_g <-> H2O2*',` | Python 看它只是字符串；字符串内部的 `*`、`<->` 由 CatMAP 自己解析。 | H2O2 吸附。 |
| 8 | `'H2O2* + * -> Ha* + OOH*',` | 字符串列表元素。 | H2O2 活化成 Ha* + OOH*。 |
| 9 | `'* + C6H12_g <-> C6H12*',` | 字符串列表元素。 | 1-己烯吸附。 |
| 10 | `'OOH* + C6H12* + * <-> OOH-C6H12* + 2* -> C6H12O* + O* + Hb*',` | 很长的字符串仍是一个列表元素，Python 不拆其中反应符号。 | 显式含 `OOH-C6H12*` 过渡态的环氧化步骤。 |
| 11 | `'C6H12O* <-> C6H12O_g + *',` | 字符串。 | 环氧产物脱附。 |
| 12 | `'Hb* + O* <-> Hb-O* + * -> OH* + *',` | 字符串。 | H 转移步骤及 `Hb-O*` 过渡态。 |
| 13 | `'Ha* + OH* -> * + H2O*',` | 字符串。 | 生成吸附水。 |
| 14 | `'H2O* <-> H2O_g + *',` | 字符串。 | 水脱附。 |
| 15 | `]` | 结束列表。 | 反应网络定义完成。 |
| 17 | `# Surfaces included in the current energies.txt data set` | 注释。 | 说明表面名称必须与 energies.txt 对应。 |
| 18 | `surface_names = [` | 列表赋值。 | 开始列出参与拟合的表面。 |
| 19 | `'tife', 'timn', 'tihf', 'tire', 'timo', 'tiv', 'tizr',` | 一行可写多个字符串列表元素。 | 第一组 Ti-M 表面。 |
| 20 | `'tico', 'titi', 'tiw', 'tita', 'ticr', 'ti', 'tini'` | 列表元素。 | 第二组表面。 |
| 21 | `]` | 结束列表。 | surface_names 完成。 |
| 23 | `# The old c_s descriptor does not exist in the present reaction network.` | 注释。 | 说明旧描述符已弃用。 |
| 24 | `# OOH_s describes H2O2 activation/intermediate stability, while O_s describes` | 注释。 | 解释 OOH_s 物理意义。 |
| 25 | `# the stability/removal of surface oxygen after epoxidation.` | 注释。 | 解释 O_s 物理意义。 |
| 26 | `descriptor_names = ['OOH_s', 'O_s']` | 创建含两个字符串的列表。 | 选择二维描述符 OOH_s 与 O_s。 |
| 28 | `# Ranges cover the present Ti-M data and leave a small margin for interpolation.` | 注释。 | 说明扫描范围选择依据。 |
| 29 | `descriptor_ranges = [[-7.0, -3.0], [-4.0, 2.5]]` | **嵌套列表**：外层列表中有两个 `[min,max]` 子列表。 | 分别规定两个描述符的能量扫描范围。 |
| 30 | `resolution = 20` | 整数赋值。 | 每个描述符方向采用 20 个网格点。 |
| 32 | `# Experimental reaction temperature: 60 degC` | 注释。 | 标明实验温度。 |
| 33 | `temperature = 333.15` | 浮点数赋值。 | 设置 333.15 K。 |
| 35 | `#` | 注释。 | 分隔。 |
| 36 | `# Species definitions / reaction conditions` | 注释。 | 物种/条件区开始。 |
| 37 | `#` | 注释。 | 分隔。 |
| 38 | `# Experimental feed:` | 注释。 | 标记实验进料说明。 |
| 39 | `#   1.2 g 1-hexene` | 注释。 | 记录 1-己烯用量。 |
| 40 | `#   1.2 g 30 wt% H2O2 aqueous solution` | 注释。 | 记录过氧化氢溶液用量。 |
| 41 | `#` | 注释。 | 分隔。 |
| 42 | `# This is a liquid-phase system...` | 注释。 | 强调 `pressure` 在这里作为有效活度参数。 |
| 43 | `# here as a first-order effective activity/reservoir parameter...` | 注释。 | 继续解释近似。 |
| 44 | `# literal gas-phase partial pressure...` | 注释。 | 说明不要把它严格理解为气相分压。 |
| 45 | `# should later include solution chemical potentials / solvation corrections.` | 注释。 | 指出未来更严格模型需加入溶剂化/液相化学势。 |
| 46 | `species_definitions = {}` | 创建空字典。 | 作为所有物种条件和元数据的容器。 |
| 48 | `# 1-hexene organic phase: pure-liquid standard-state approximation, a ~ 1` | 注释。 | 说明己烯活度近似。 |
| 49 | `species_definitions['C6H12_g'] = {` | `dict[key] = value` 可逐步往已有字典里添加项目；value 还能是另一个字典。 | 开始定义 C6H12_g 条件。 |
| 50 | `'pressure': 1.0` | 内层字典键值对。 | 有效活度设为 1。 |
| 51 | `}` | 结束内层字典。 | C6H12_g 定义完成。 |
| 53 | `# 30 wt% H2O2 aqueous solution corresponds approximately to` | 注释。 | 说明后面的摩尔分数来源。 |
| 54 | `# x(H2O2) = 0.185 and x(H2O) = 0.815.` | 注释。 | 给出 H2O2/H2O 有效比例。 |
| 55 | `species_definitions['H2O2_g'] = {` | 往字典增加一个嵌套字典。 | 定义 H2O2_g。 |
| 56 | `'pressure': 0.185` | 内层键值对。 | H2O2 有效活度 0.185。 |
| 57 | `}` | 结束内层字典。 | H2O2_g 完成。 |
| 59 | `species_definitions['H2O_g'] = {` | 增加字典项。 | 定义 H2O_g。 |
| 60 | `'pressure': 0.815` | 键值对。 | H2O 有效活度 0.815。 |
| 61 | `}` | 结束字典。 | H2O_g 完成。 |
| 63 | `# Epoxide is absent at the beginning of the experiment...` | 注释。 | 说明初始产物近似为 0。 |
| 64 | `# value is used instead of exactly zero for numerical robustness.` | 注释。 | 说明为什么不用严格 0。 |
| 65 | `species_definitions['C6H12O_g'] = {` | 增加字典项。 | 定义产物储库。 |
| 66 | `'pressure': 1e-20` | `1e-20` 是科学计数法浮点数。 | 用极小正数近似“没有初始产物”。 |
| 67 | `}` | 结束字典。 | 产物定义完成。 |
| 69 | `# Surface site` | 注释。 | 开始位点定义。 |
| 70 | `species_definitions['s'] = {` | 嵌套字典。 | 定义抽象位点 `s`。 |
| 71 | `'site_names': ['111'],` | value 是一个单元素列表。 | 把 s 映射到 111 位点。 |
| 72 | `'total': 1` | 整数键值对。 | 总位点数归一化为 1。 |
| 73 | `}` | 结束字典。 | 位点定义完成。 |
| 75 | `# Ha/Hb denote H atoms...` | 注释。 | 说明 Ha/Hb 非标准化学式。 |
| 76 | `# chemical formulae, so their elemental compositions must be supplied explicitly.` | 注释。 | 说明需要手工提供组成。 |
| 77 | `species_definitions['Ha_s'] = {` | 增加嵌套字典。 | 定义 Ha_s。 |
| 78 | `'composition': {'H': 1}` | 内层又嵌字典。 | 指定 Ha_s 含 1 个 H。 |
| 79 | `}` | 结束字典。 | Ha_s 完成。 |
| 81 | `species_definitions['Hb_s'] = {` | 增加字典项。 | 定义 Hb_s。 |
| 82 | `'composition': {'H': 1}` | 嵌套字典。 | Hb_s = H1。 |
| 83 | `}` | 结束字典。 | Hb_s 完成。 |
| 85 | `species_definitions['Hb-O_s'] = {` | 增加字典项。 | 定义 Hb-O_s。 |
| 86 | `'composition': {'H': 1, 'O': 1}` | 内层字典含两个元素计数。 | Hb-O_s = H1O1。 |
| 87 | `}` | 结束字典。 | Hb-O_s 完成。 |
| 89 | `#` | 注释。 | 分隔。 |
| 90 | `# Data / parser parameters` | 注释。 | 数据与解析参数区。 |
| 91 | `#` | 注释。 | 分隔。 |
| 92 | `data_file = 'alkene_epoxidation.pkl'` | 字符串赋值；`.pkl` 常表示 Python 序列化文件。 | 指定 CatMAP 数据缓存文件。 |
| 93 | `input_file = 'energies.txt'` | 字符串赋值。 | 指定真正的多表面能量输入表。 |
| 95 | `#` | 注释。 | 分隔。 |
| 96 | `# Thermodynamic parameters` | 注释。 | 热力学配置区。 |
| 97 | `#` | 注释。 | 分隔。 |
| 98 | `# The experiment is liquid phase...` | 注释。 | 说明 `_g` 只是储库标签。 |
| 99 | `# Until explicit liquid/solvation free energies are available...` | 注释。 | 说明当前近似原因。 |
| 100 | `# adding a gas-phase entropy correction...` | 注释。 | 避免不恰当气相熵修正。 |
| 101 | `gas_thermo_mode = 'frozen_gas'` | 字符串配置赋值。 | 储库物种采用 frozen_gas。 |
| 102 | `adsorbate_thermo_mode = 'frozen_adsorbate'` | 字符串配置赋值。 | 吸附物采用 frozen_adsorbate。 |
| 103 | `pressure_mode = 'static'` | 字符串配置赋值。 | 边界活度保持静态。 |
| 105 | `#` | 注释。 | 分隔。 |
| 106 | `# Scaling parameters` | 注释。 | 标度关系参数区。 |
| 107 | `#` | 注释。 | 分隔。 |
| 108 | `# Force the two selected descriptors...` | 注释。 | 解释约束目的。 |
| 109 | `# adsorption energies. Other species are fitted...` | 注释。 | 其他物种交给广义线性 scaler。 |
| 110 | `scaling_constraint_dict = {` | 创建字典。 | 开始定义 scaling 约束。 |
| 111 | `'OOH_s': ['+', 0, None],` | value 是列表，可混合字符串、整数和 `None`。 | 约束 OOH_s 主要对应第一描述符。 |
| 112 | `'O_s': [0, '+', None],` | 同上。 | 约束 O_s 主要对应第二描述符。 |
| 113 | `}` | 结束字典。 | scaling 约束完成。 |
| 115 | `#` | 注释。 | 分隔。 |
| 116 | `# Solver parameters` | 注释。 | 求解器参数区。 |
| 117 | `#` | 注释。 | 分隔。 |
| 118 | `decimal_precision = 100` | 整数赋值。 | 高精度运算使用 100 位精度。 |
| 119 | `tolerance = 1e-50` | 科学计数法浮点数。 | 设置极严格收敛容差。 |
| 120 | `max_rootfinding_iterations = 100` | 整数赋值。 | 根求解最多 100 次迭代。 |
| 121 | `max_bisections = 3` | 整数赋值。 | 最多 3 次二分/细分辅助。 |

---

# 5. `2-creating_microkinetic_model/test.ipynb`

## Cell 1：运行模型并绘图

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `import page` | 普通模块导入；如果模块不存在会报 `ModuleNotFoundError`。 | 当前后续代码没有使用 `page`，很可能是旧项目残留。 |
| 3 | `# -*- coding: utf-8 -*-` | 这是历史上常见的源文件编码声明；在现代 Python 3 的 UTF-8 环境里通常不必写。 | 对当前 Notebook 基本没有实际作用。 |
| 4 | `# page.encoding='utf-8'` | 被 `#` 注释掉的代码不会执行。 | 旧设置，当前无效。 |
| 6 | `from catmap import ReactionModel` | 指定对象导入。 | 导入 CatMAP 模型类。 |
| 8 | `mkm_file = 'alkene_epoxidation.mkm'` | 字符串赋值。 | 保存模型配置文件名。 |
| 9 | `model = ReactionModel(setup_file=mkm_file)` | 类实例化；`setup_file=...` 是关键字参数。 | 读取 `.mkm` 并构造完整模型。 |
| 10 | `model.output_variables += ['production_rate','rate','rate_control','coverage','selectivity_control','rxn_order']` | `+=` 对列表相当于扩展列表；右侧也是列表。 | 追加需要 CatMAP 计算和保存的输出变量。 |
| 11 | `model.run()` | 调用对象方法。 | 在描述符空间真正求解稳态模型。 |
| 13 | `from catmap import analyze` | 导入模块/对象。 | 获取 CatMAP 后处理模块。 |
| 14 | `vm = analyze.VectorMap(model)` | 通过模块中的类构造对象。 | 创建二维描述符绘图对象。 |
| 15 | `vm.plot_variable = 'rate'` | 对象属性赋值。 | 指定绘制基元反应速率。 |
| 16 | `vm.log_scale = True` | `True` 是布尔值。 | 速率图采用对数尺度。 |
| 17 | `vm.min = 1e-25` | 科学计数法赋值。 | 设置绘图下限。 |
| 18 | `vm.max = 1e3` | 科学计数法赋值。 | 设置绘图上限。 |
| 19 | `vm.plot(save='rate.pdf')` | 方法调用 + 关键字参数。 | 保存 `rate.pdf`。 |
| 21 | `vm.unique_only = False` | `False` 是布尔值。 | 暂时允许显示所有反应项。 |
| 22 | `vm.plot(save='all_rates.pdf')` | 方法调用。 | 保存全部速率图。 |
| 23 | `vm.unique_only = True` | 布尔属性赋值。 | 恢复只保留唯一项。 |
| 25 | `vm.production_rate_map = model.production_rate_map` | 可以把一个对象的属性直接赋给另一个对象的属性。 | 把模型已有生产速率图数据挂给绘图对象。 |
| 26 | `vm.threshold = 1e-30` | 浮点属性赋值。 | 设置绘图阈值。 |
| 27 | `vm.plot_variable = 'production_rate'` | 字符串属性赋值。 | 切换到生产速率。 |
| 28 | `vm.plot(save='production_rate.pdf')` | 方法调用。 | 保存生产速率图。 |
| 30 | `vm.descriptor_labels = ['H reactivity [eV]', 'C2H5 reactivity [eV]']` | 给属性赋字符串列表。 | **与当前 `OOH_s/O_s` 描述符不一致，是旧示例残留。** |
| 31 | `vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}` | 字典常用来集中保存“关键字参数”；`kwargs` 通常就是 keyword arguments 的缩写。 | 设置 Matplotlib 子图边距。 |
| 32 | `vm.plot(save='pretty_production_rate.pdf')` | 方法调用。 | 保存美化后的生产速率图。 |
| 34 | `vm.plot_variable = 'coverage'` | 属性赋值。 | 切换到覆盖度。 |
| 35 | `vm.log_scale = False` | 布尔值。 | 覆盖度不用对数轴。 |
| 36 | `vm.min = 0` | 整数属性赋值。 | 覆盖度下限 0。 |
| 37 | `vm.max = 1` | 整数属性赋值。 | 覆盖度上限 1。 |
| 38 | `vm.plot(save='coverage.pdf')` | 方法调用。 | 保存覆盖度图。 |
| 40 | `vm.include_labels = ['C2H5_s']` | 列表属性赋值。 | **当前模型没有 C2H5_s，属于旧网络残留。** |
| 41 | `vm.plot(save='C2H5_coverage.pdf')` | 方法调用。 | 试图绘制 C2H5_s 覆盖度，因此当前代码逻辑需要更新。 |
| 43 | `sa = analyze.ScalingAnalysis(model)` | 类实例化。 | 创建标度关系分析对象。 |
| 44 | `sa.plot(save='scaling.pdf')` | 方法调用。 | 保存 scaling 图。 |

### Cell 1 建议

把：

```python
vm.descriptor_labels = ['H reactivity [eV]', 'C2H5 reactivity [eV]']
vm.include_labels = ['C2H5_s']
```

更新为与当前机理一致的 `OOH_s / O_s` 相关标签。

## Cell 2：从 `.log` 恢复模型并导出表格

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `from ase.calculators.calculator import InputError` | 从深层子模块导入异常类。 | 用一个明确异常表示输入日志文件不唯一。 |
| 2 | `from glob import glob` | `glob` 是标准库的通配符文件搜索函数。 | 搜索当前目录的 `*.log`。 |
| 3 | `import sys` | 导入标准库模块。 | 当前单元没用到 `sys`，可删除。 |
| 4 | `from catmap.model import ReactionModel` | 导入类。 | 用日志恢复模型。 |
| 6 | `model.output_variables += [...]` | 列表增强赋值。 | 追加输出变量。**但若独立执行该 Cell，此时 `model` 尚未定义，会 `NameError`。** |
| 8 | `output_variable = 'production_rate'` | 字符串赋值。 | 指定导出生产速率。 |
| 9 | `logfile = glob('*.log')` | 函数调用返回列表。 | 获取所有日志文件路径。 |
| 10 | `if len(logfile) > 1:` | `len()` 求长度；`>` 比较大小。 | 多于一个日志时认为有歧义。 |
| 11 | `raise InputError('Ambiguous logfile...')` | 主动抛异常。 | 阻止随机选择错误日志。 |
| 12 | `model = ReactionModel(setup_file=logfile[0])` | 列表索引从 0 开始；`logfile[0]` 是第一个元素。 | 从第一个 log 恢复模型。若没有 log，这里会 `IndexError`。 |
| 14 | `if output_variable == 'rate_control':` | 字符串相等比较。 | rate_control 是二维输出，需要特殊处理。 |
| 15 | `dim = 2` | 整数赋值。 | 标记二维。 |
| 16 | `else:` | `else` 在 `if` 不成立时执行。 | 非 rate_control 走普通分支。 |
| 17 | `dim = 1` | 整数赋值。 | 标记一维。 |
| 19 | `labels = model.output_labels['production_rate']` | 对象属性中取字典，再按 key 取值。 | 获取生产速率列标签；这里写死了 key。 |
| 21 | `def flatten_2d(output):` | 定义辅助函数。 | 把二维结果展开成一维。 |
| 22 | `"Helper function for flattening rate_control output"` | 这是单引号/双引号字符串作为函数首句，同样可作为 docstring。 | 说明函数用途。 |
| 23 | `flat = []` | 空列表。 | 保存展开结果。 |
| 24 | `for x in output:` | 遍历外层列表。 | 每次取一行二维结果。 |
| 25 | `flat += x` | 列表 `+=` 会把右侧可迭代对象的元素扩展进左侧列表。 | 把一行元素接到 flat 末尾。 |
| 26 | `return flat` | 返回值。 | 返回一维列表。 |
| 28 | `#flatten rate_control labels` | 注释。 | 下面处理二维标签。 |
| 29 | `if output_variable == 'rate_control':` | 条件判断。 | 只对 rate_control 执行。 |
| 30 | `flat_labels = []` | 空列表。 | 保存新标签。 |
| 31 | `for i in labels[0]:` | 列表索引 + 循环。 | 遍历第一维标签。 |
| 32 | `for j in labels[1]:` | 嵌套循环。 | 对每个 i 再遍历第二维 j。 |
| 33 | `flat_labels.append('d'+i+'/d'+j)` | `+` 可拼接字符串；`append()` 加入列表。 | 生成 `dA/dB` 形式控制系数名。 |
| 34 | `labels = flat_labels` | 变量重新绑定。 | 用扁平标签替换原标签。 |
| 36 | `#flatten elementary-step specific labels` | 注释。 | 下面处理基元步骤标签。 |
| 37 | `if output_variable in ['rate','rate_constant','forward_rate_constant','reverse_rate_constant']:` | `in` 可检查某值是否存在于列表中。 | 只有反应步骤类输出需要把结构化标签转成反应字符串。 |
| 38 | `str_labels = []` | 空列表。 | 保存字符串标签。 |
| 39 | `for label in labels:` | 循环。 | 逐个处理反应标签。 |
| 40 | `states = ['+'.join(s) for s in label]` | **列表推导式** `[表达式 for 变量 in 可迭代对象]` 是快速生成列表的语法。 | 把每个反应状态中的物种用 `+` 拼接。 |
| 41 | `if len(states) == 2:` | 长度判断。 | 两态反应单独格式化。 |
| 42 | `new_label = '<->'.join(states)` | `join()` 连接字符串列表。 | 生成 `A<->B` 标签。 |
| 43 | `else:` | 否则分支。 | 处理含显式 TS 的三态。 |
| 44 | `new_label = states[0]+'<->'+states[1]+'->'+states[2]` | 字符串拼接 + 列表索引。 | 生成 `初态<->TS->终态`。 |
| 45 | `str_labels.append(new_label)` | 列表追加。 | 保存新标签。 |
| 46 | `labels = str_labels` | 重新赋值。 | 替换原结构化标签。 |
| 48 | `table = '\t'.join(list(['descriptor-'+d for d in model.descriptor_names])+list(labels))+'\n'` | 一行组合了列表推导式、类型转换 `list()`、列表拼接、`join()` 和字符串拼接。 | 构造导出表头：描述符列 + 输出列。 |
| 50 | `for pt, output in getattr(model,output_variable+'_map'):` | `getattr(obj, name)` 按**字符串动态访问属性**；这里属性名由 `output_variable` 拼出来。 | 动态读取如 `production_rate_map`，避免写死多个 if。 |
| 51 | `if dim == 2:` | 条件判断。 | 二维输出需要摊平。 |
| 52 | `output = flatten_2d(output)` | 自定义函数调用。 | 把二维矩阵转一维。 |
| 53 | `table += '\t'.join([str(float(i)) for i in pt+output])+'\n'` | 列表推导式 + `float()` + `str()` + `join()` + `+=`。 | 把描述符坐标与结果转成一行文本并追加。 |
| 55 | `f = open(output_variable+'_table.txt','w')` | 直接 `open()` 会返回文件对象；不用 `with` 时必须手工关闭。 | 创建如 `production_rate_table.txt`。 |
| 56 | `f.write(table)` | 文件写入。 | 写整个表格。 |
| 57 | `f.close()` | 显式关闭文件。 | 刷新缓冲并释放文件资源。 |

### Cell 2 当前三个明显问题

1. `model.output_variables += ...` 出现在 `model = ReactionModel(...)` **之前**，独立运行会 `NameError`。
2. 只检查 `len(logfile) > 1`，没检查 `len(logfile) == 0`。
3. `labels = model.output_labels['production_rate']` 写死了 key，应考虑改成 `model.output_labels[output_variable]`。

## Cell 3：读取并打印覆盖度

| 行号 | 代码 | Python 基础概念（先学这个） | 本项目中的作用 |
|---:|---|---|---|
| 1 | `from catmap.model import ReactionModel` | 导入类。 | 获取模型恢复功能。 |
| 3 | `model = ReactionModel(setup_file='alkene_epoxidation.log')` | 类实例化 + 关键字参数。 | 从已运行的 log 恢复模型。 |
| 5 | `#for MgO, cvgs in model.coverage_map:` | 被注释掉的代码不会执行。 | 旧示例循环，保留作参考。 |
| 6 | `# print('descriptors:', MgO)` | 注释掉的函数调用。 | 旧调试输出。 |
| 7 | `#print('coverages', cvgs)` | 注释掉的函数调用。 | 旧覆盖度输出。 |
| 9 | `labels = model.output_labels['coverage']` | 对象属性 + 字典索引。 | 获取覆盖度各分量对应的物种名。 |
| 10 | `for MgO, cvg in model.coverage_map:` | `for` 循环 + 两变量解包。 | 遍历“描述符坐标 + 覆盖度向量”。变量名 `MgO` 与当前 Ti-M 项目无关。 |
| 11 | `print('descriptors', MgO)` | 输出函数。 | 打印描述符坐标。 |
| 12 | `print('intermediates', labels)` | 输出函数。 | 打印物种标签。 |
| 13 | `print('coverages', [float(c) for c in cvg])` | 列表推导式 + 类型转换。 | 把高精度数值转成普通 float 后打印覆盖度。 |

---

# 6. 把这些 Python 基础串起来

这个项目实际上反复使用了十几个核心 Python 概念：

| 基础概念 | 最典型代码 | 你需要掌握的核心点 |
|---|---|---|
| 变量赋值 | `E0 = energy` | 名字绑定到对象；后续可以重新赋值。 |
| 字典 | `ref_dict = {'O': -4.37}` | 最适合保存“名称 → 参数/数据”的映射。 |
| 列表 | `surface_names = [...]` | 有序可变容器。 |
| 元组 | `('s',)` | 有序不可变容器；单元素元组必须有逗号。 |
| 条件 | `if site == '111':` | 只在条件为真时执行缩进块。 |
| 循环 | `for atom in composition:` | 对一系列元素重复执行逻辑。 |
| 函数 | `def get_formation_energies(...):` | 把重复逻辑封装为可复用单元。 |
| 返回值 | `return formation_energies` | 函数把结果交给调用者。 |
| 字符串方法 | `split / replace / join` | 本项目命名解析和表格生成的核心。 |
| 异常 | `try / except / raise` | 让错误显式暴露，而不是悄悄产生错误结果。 |
| 文件 IO | `with open(...)` | 读取/写入磁盘文件。 |
| 对象属性 | `rxm.surface_names = ...` | CatMAP 大量配置都通过对象属性完成。 |
| 列表推导式 | `[float(c) for c in cvg]` | 简洁地“遍历 + 转换 + 生成新列表”。 |
| 动态属性 | `getattr(model, name)` | 当属性名由字符串动态决定时使用。 |

---

# 7. 三个核心数据对象

## `abinitio_energies`

最原始的总能数据库：

```text
key = 物种名 + site
value = DFT 总能 / 自由能
```

## `formation_energies`

由：

```python
formation_energies = get_formation_energies(abinitio_energies, ref_dict)
```

生成。

当前定义可概括为：

```text
气相：
E_form = E_gas - Σ n_i E_ref,i

表面态 / 过渡态：
E_form = E_system - E_slab - Σ n_i E_ref,i
```

## `species_definitions`

CatMAP 模型层的物种元数据容器，主要保存：

- 储库有效活度/pressure；
- 位点映射；
- 特殊物种元素组成；
- 其他模型属性。

---

# 8. 当前代码中最值得优先理解/修改的地方

1. `surface = None if site == 'gas' else 'tini'` 把表面名**硬编码为 tini**；若以后批量处理 TiFe、TiMn、TiHf 等，应把 `surface_name` 变成函数参数。
2. `test.ipynb` 中 `H reactivity`、`C2H5 reactivity`、`C2H5_s` 与当前 `OOH_s / O_s` 描述符网络不一致，是旧代码残留。
3. `import page`、变量名 `MgO` 也属于旧项目痕迹。
4. Notebook 的 Cell 2 对 `model` 的使用顺序有问题，且缺少“没有 `.log` 文件”的分支。

如果后面继续学习，建议顺序是：**字典 → 循环/条件 → 函数 → 字符串处理 → 异常 → 文件 IO → 类与对象属性 → 列表推导式/动态属性**。这些已经覆盖了当前项目绝大多数 Python 语法。
