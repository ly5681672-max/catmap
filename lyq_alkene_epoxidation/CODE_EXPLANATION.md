# `lyq_alkene_epoxidation` 代码伴读手册
## Python 基础 + CatMAP 微观动力学 + 变量追踪版

> 本文不再采用“每一行都重复解释一次”的阅读方式，而采用两层结构：
>
> **文件级：固定 5 步模板**  
> ① 一句话总览 → ② 整体架构/流程 → ③ Python 基础概念 → ④ 关键变量&核心逻辑 → ⑤ 记忆要点+易错坑
>
> **代码块级：D 双轨科研伴读法**  
> 当前任务 → 参数卡 → 小段代码 → Python 在做什么 → 科研/计算上在做什么 → 真实变量追踪 → 状态快照
>
> 阅读目标不是“背住所有变量”，而是做到：**任何时候看到一个变量，都能马上知道它从哪里来、现在是什么、接下来去哪里。**

---

# 0. 先建立整个项目的认知地图

## 【1. 一句话总览】

这个项目把 **DFT/自由能数据** 转换成 CatMAP 可读取的形成能表，再建立 **1-己烯 + H2O2 环氧化微观动力学模型**，最终计算速率、覆盖度、生产速率、控制系数等结果。

---

## 【2. 整体架构 / 数据流】

```text
DFT / 自由能原始数据
        │
        ▼
data_only_not_executable.py
generate_input1.py
        │
        │ 形成能换算
        ▼
*_energies.txt / *_energies.csv
        │
        │ 汇总不同 Ti-M 表面
        ▼
energies.txt
        │
        ▼
alkene_epoxidation.mkm
        │
        │ 反应网络 + 描述符 + 条件 + 求解参数
        ▼
CatMAP ReactionModel
        │
        ▼
model.run()
        │
        ├── rate
        ├── production_rate
        ├── coverage
        ├── rate_control
        ├── selectivity_control
        └── rxn_order
        │
        ▼
test.ipynb
        │
        ├── 二维图
        ├── 输出表格
        └── 覆盖度检查
```

你阅读任何代码时，先判断自己目前处于哪一层：

```text
A. 原始能量
B. 形成能
C. CatMAP 输入表
D. 微观动力学配置
E. 模型求解
F. 后处理/绘图
```

---

## 【3. 本项目核心 Python 基础概念】

| Python 概念 | 典型代码 | 在本项目中的作用 |
|---|---|---|
| 变量赋值 | `E0 = energy` | 保存当前计算状态 |
| 字典 `dict` | `ref_dict = {...}` | 保存“名称 → 数值/参数”映射 |
| 列表 `list` | `surface_names = [...]` | 保存有序的一组表面/物种 |
| 元组 `tuple` | `('s',)` | CatMAP 中保存固定名称集合 |
| `for` 循环 | `for atom in composition:` | 逐个物种/原子执行相同计算 |
| `if` 条件 | `if site == '111':` | 区分气相与表面态 |
| 函数 `def` | `get_formation_energies()` | 封装形成能计算 |
| `return` | `return formation_energies` | 把函数结果交给外部 |
| 字符串处理 | `split/replace/join` | 解析物种名、生成表格 |
| 异常 | `try/except/raise` | 让命名或解析错误显式暴露 |
| 文件 IO | `with open(...)` | 写 TXT / CSV |
| 类与对象 | `ReactionModel()` | 构造 CatMAP 模型 |
| 对象属性 | `model.temperature = ...` | 配置/读取模型状态 |
| 列表推导式 | `[float(c) for c in cvg]` | 批量转换数据 |
| 动态属性 | `getattr(model, name)` | 根据字符串决定读取哪个结果 |

### Python 学习规则

第一次遇到某个知识点：**详细讲。**

之后再次出现：

> ♻ **Python 复习：** 已经见过的“字典遍历 + 两变量解包”。

这样保留复习效果，但不让重复解释淹没计算主线。

---

## 【4. 全项目变量速查表】

这是本文最重要的一张表。**忘记变量时不要往前翻几十行，直接查这里。**

| 变量 | 类型 | 一句话含义 | 从哪里来 | 去哪里 |
|---|---|---|---|---|
| `abinitio_energies` | `dict` | 原始 DFT/自由能数据库 | 手工输入 | `get_formation_energies()` |
| `ref_dict` | `dict` | H/C/O 与裸板参考能 | 手工输入 + 裸板能 | 形成能计算 |
| `formula_map` | `dict` | 非标准名字 → 真正化学式 | 手工定义 | ASE 化学式解析 |
| `energy_dict` | 函数参数 | 传入函数的能量字典 | 实参 `abinitio_energies` | 函数内部循环 |
| `references` | 函数参数 | 传入函数的参考能字典 | 实参 `ref_dict` | 扣参考能 |
| `key` | `str` | 当前字典键，如 `OOH_111` | `dict.items()` | 拆成 `name/site` |
| `energy` | `float` | 当前物种原始能量 | `dict.items()` | 赋给 `E0` |
| `name` | `str` | 当前物种名，如 `OOH` | `key.split()` | 化学式解析/输出 |
| `site` | `str` | 当前位点，如 `111/gas` | `key.split()` | 判断是否减裸板 |
| `E0` | `float` | 当前正在修改的“工作能量” | `energy` | 最终形成能 |
| `formula` | `str` | ASE 可识别化学式 | `formula_map`/`replace()` | `string2symbols()` |
| `composition` | `list` | 元素序列，如 `['O','O','H']` | ASE | 逐原子扣参考能 |
| `formation_energies` | `dict` | 最终形成能数据库 | 函数计算 | TXT/CSV |
| `frequency_dict` | `dict` | 各物种振动频率 | 当前多为空 | TXT/CSV |
| `file_name` | `str` | 输出文件名 | 手工设置 | `open()` / parser |
| `rxm` | 对象 | 最小 CatMAP 模型 | `ReactionModel()` | `TableParser` |
| `parser` | 对象 | CatMAP 表格解析器 | `TableParser(rxm)` | `parse()` |
| `rxn_expressions` | `list` | 基元反应网络 | `.mkm` | CatMAP |
| `surface_names` | `list` | 参与建模的 Ti-M 表面 | `.mkm` | scaler/model |
| `descriptor_names` | `list` | 二维描述符名称 | `.mkm` | 描述符空间 |
| `descriptor_ranges` | `list` | 描述符扫描范围 | `.mkm` | 网格 |
| `species_definitions` | `dict` | 活度、位点、元素组成 | `.mkm` | CatMAP |
| `data_file` | `str` | CatMAP 缓存文件 | `.mkm` | 模型保存/恢复 |
| `input_file` | `str` | 多表面能量文件 | `.mkm` | TableParser |
| `model` | 对象 | 完整微观动力学模型 | `ReactionModel(setup_file=...)` | `run()`/分析 |
| `vm` | 对象 | 二维结果绘图器 | `VectorMap(model)` | PDF 图 |
| `sa` | 对象 | 标度关系分析器 | `ScalingAnalysis(model)` | scaling 图 |
| `labels` | list-like | 输出各列对应名称 | `model.output_labels[...]` | 导出/打印 |
| `cvg` | list-like | 某描述符点的覆盖度 | `coverage_map` | 打印/后处理 |

### 变量名快速判断口诀

```text
xxx_dict        → 映射表
xxx_names       → 名称列表
xxx_ranges      → 范围
xxx_mode        → 计算模式
xxx_file        → 文件路径/文件名
xxx_map         → 某个网格上的结果
rxm / model     → CatMAP 模型对象
parser          → 解析输入文件
vm              → VectorMap，负责画二维图
E0              → 当前正在加工的能量工作变量
```

---

## 【5. 总体记忆要点 + 易错坑】

### 记忆口诀

> **总能进来 → 减参考 → 形成能 → 写表 → 定机理 → 跑 CatMAP → 读 map。**

### 当前最容易混淆的 4 组变量

```text
abinitio_energies  ≠ formation_energies
原始总能              参考态换算后的形成能

energy ≠ E0
只读原始值     当前不断被扣参考能的工作值

name ≠ formula
项目物种名     ASE 真正解析的化学式

rxm ≠ model
输入解析阶段模型   正式微观动力学模型
```

---

# 1. `1-generating_input_file/data_only_not_executable.py`

# 文件级固定 5 步

## 【1. 一句话总览】

把一个表面的原始 DFT/自由能转换为 **CatMAP 形成能 TXT**，随后创建最小 `ReactionModel` 验证该 TXT 能否被 CatMAP 正确解析。

## 【2. 整体架构 / 流程】

```text
abinitio_energies
        │
        ├── ref_dict
        └── formula_map
        │
        ▼
get_formation_energies()
        │
        ▼
formation_energies
        │
        ├── print 检查
        └── frequency_dict
        │
        ▼
make_input_file()
        │
        ▼
tini_energies.txt
        │
        ▼
ReactionModel + TableParser
        │
        ▼
parser.parse()
```

## 【3. 核心 Python 基础】

本文件重点复习：

- `dict`：能量数据库/参考能数据库/名称映射；
- `for` + `.items()`：遍历物种；
- `split()`：从 `OOH_111` 中拆出物种和位点；
- `dict.get()`：有特殊映射用映射，没有就用默认值；
- `try/except/raise`：化学式解析失败时停止；
- `with open()`：写入 TXT；
- 类实例化与对象属性：`ReactionModel()`、`rxm.surface_names = ...`。

## 【4. 关键变量 & 核心逻辑】

真正核心只有两个算法：

```text
算法 A：形成能换算
算法 B：把形成能组织成 CatMAP TableParser 所需的表格
```

形成能当前定义：

```text
气相：
E_form = E_gas - Σ n_i E_ref,i

表面态 / TS：
E_form = E_system - E_slab - Σ n_i E_ref,i
```

## 【5. 记忆要点 + 易错坑】

> **先拆名字，再判 gas/111，再减参考。**

最关键的坑：

```python
name, site = key.split('_', 1)
```

这里必须保留 `1`。如果以后物种名中还有下划线，直接 `split('_')` 可能产生超过两个结果，引发：

```text
ValueError: too many values to unpack
```

---

# 1.1 输入数据库：程序最开始“知道了什么”

## 当前任务

> 🎯 **把所有原始能量和参考态装进内存。这里还没有开始计算形成能。**

### 参数卡

```text
abinitio_energies
├─ 类型：dict[str, float]
├─ key：物种名_site
├─ value：原始 DFT/自由能
└─ 单位：eV

ref_dict
├─ 类型：dict[str, float]
├─ H/C/O：原子参考能
└─ 111：裸板参考能

formula_map
├─ 类型：dict[str, str]
└─ 解决 Ha/Hb 等不是正规元素符号的问题
```

### 代码

```python
from ase.symbols import string2symbols

abinitio_energies = {
    'H2O2_gas': -18.122091,
    ...
    'OOH_111': -877.1433925,
    ...
    'slab_111': -863.25438,
}

ref_dict = {
    'O': -4.37,
    'C': -9.28,
    'H': -1.11,
    '111': abinitio_energies['slab_111'],
}

formula_map = {
    'Ha': 'H',
    'Hb': 'H',
    'Hb-O': 'HO',
}
```

### 🐍 Python 在做什么？

`{key: value}` 创建字典。

```python
abinitio_energies['slab_111']
```

表示按 key 查询字典中的值。

```python
from ase.symbols import string2symbols
```

这行代码不只是“从模块里导入一个函数”，还需要知道**具体导入的函数是什么、它能做什么**。

### ① `ase.symbols` 是什么？

`ase` 是 **Atomic Simulation Environment（原子模拟环境）** 这个 Python 科学计算库；`ase.symbols` 是 ASE 中专门处理**元素符号和化学式字符串**的模块。

### ② `string2symbols` 是什么函数？

`string2symbols()` 是 `ase.symbols` 模块提供的一个**化学式解析函数**。

它的任务是：

> **把一个化学式字符串展开成“逐个原子的元素符号列表”。**

例如：

```python
string2symbols('H2O2')
```

返回：

```python
['H', 'H', 'O', 'O']
```

再例如本项目中的：

```python
string2symbols('OOH')
```

返回：

```python
['O', 'O', 'H']
```

因此它的输入和输出可以记成：

```text
输入：化学式字符串 str
      例如 'H2O2'、'OOH'、'C6H12O'

             ↓ string2symbols()

输出：元素符号列表 list[str]
      例如 ['H','H','O','O']
```

### ③ 为什么本项目必须用这个函数？

因为后面形成能计算需要**一个原子一个原子地减去参考能**：

```python
for atom in composition:
    E0 -= references[atom]
```

所以程序不能只知道：

```text
formula = 'OOH'
```

还必须把它展开成：

```text
composition = ['O', 'O', 'H']
```

这样循环才能依次执行：

```text
第1次：atom = O → 减一次 O 参考能
第2次：atom = O → 再减一次 O 参考能
第3次：atom = H → 减一次 H 参考能
```

也就是说，这个导入函数在本项目中的角色是：

```text
化学式字符串
    ↓
string2symbols()
    ↓
逐原子元素列表
    ↓
for atom in composition
    ↓
逐原子扣除参考能
```

### ④ `from ... import ...` 在这里具体意味着什么？

```python
from ase.symbols import string2symbols
```

可以拆成：

```text
ase.symbols
= 函数所在的模块

string2symbols
= 从这个模块中取出来的具体函数
```

由于已经把 `string2symbols` 这个函数直接导入当前文件，所以后面可以直接写：

```python
string2symbols(formula)
```

而不用写完整路径：

```python
ase.symbols.string2symbols(formula)
```

### ⑤ 一句话记忆

> **`string2symbols` = “化学式字符串 → 一个个元素符号”，为后面的逐原子参考能扣除做准备。**

### 🧪 科研上在做什么？

这里定义了三套完全不同的数据：

```text
abinitio_energies
= 实际体系计算出来的总能

ref_dict
= 形成能零点

formula_map
= 名称翻译器
```

### 状态快照

```text
此时已经有：
✓ 原始体系能量
✓ 元素参考能
✓ 裸板参考能
✓ 特殊名字映射

此时还没有：
× formation energy
× CatMAP 输入文件
× 微观动力学结果
```

---

# 1.2 `get_formation_energies()`：本文件最重要的函数

## 函数卡

| 项目 | 内容 |
|---|---|
| 函数 | `get_formation_energies(energy_dict, references)` |
| 一句话作用 | 原始总能 → 形成能 |
| 输入 1 | `energy_dict`，实际调用时就是 `abinitio_energies` |
| 输入 2 | `references`，实际调用时就是 `ref_dict` |
| 输出 | `formation_energies` |
| 是否修改输入 | 否 |
| 核心 Python | 函数、字典循环、字符串拆分、异常、返回值 |
| 核心科研逻辑 | 扣裸板 + 扣元素参考能 |

### 形参和实参不要混

定义函数时：

```python
def get_formation_energies(energy_dict, references):
```

这里的：

```text
energy_dict
references
```

叫 **形参**。

真正调用：

```python
formation_energies = get_formation_energies(
    abinitio_energies,
    ref_dict
)
```

这时候对应关系是：

```text
energy_dict  ← abinitio_energies
references   ← ref_dict
```

以后只要看到 `energy_dict`，脑中直接替换成：

> “传进来的那一整套原始能量。”

---

## 步骤 A：创建结果容器

```python
formation_energies = {}
```

### 🐍 Python

`{}` 是空字典。

### 🧪 当前计算

只是准备一个空容器：

```text
formation_energies = {}
```

后面每算完一个物种，就加入一项。

---

## 步骤 B：一次取出一个物种

```python
for key, energy in energy_dict.items():
```

### 🐍 Python

`.items()` 每次给出一对：

```text
(key, value)
```

`key, energy` 是**序列解包**。

### 🧪 真实变量追踪：以 OOH 为例

当循环走到：

```python
'OOH_111': -877.1433925
```

此刻：

```text
key    = 'OOH_111'
energy = -877.1433925 eV
```

### 状态快照

```text
现在电脑只选中了“一个”物种：
OOH_111

下一步：
把 OOH 和 111 拆开
```

---

## 步骤 C：拆物种名与位点

```python
if '_' not in key:
    raise ValueError(...)

name, site = key.split('_', 1)
```

### 🐍 Python

```python
split('_', 1)
```

意思是：

> 以 `_` 为分隔符，**最多只切 1 次**。

所以：

```text
OOH_111
  │
  ▼
['OOH', '111']
  │       │
 name    site
```

### 🧪 当前计算

```text
key  = OOH_111
name = OOH
site = 111
```

这两个变量以后分工不同：

```text
name → 判断是什么化学物种
site → 判断是不是表面态
```

---

## 步骤 D：跳过裸板

```python
if 'slab' in name:
    continue
```

### 🐍 Python

`continue`：

> 结束当前这一轮循环，直接处理下一个物种。

### 🧪 科研逻辑

裸板：

```text
slab_111
```

只用于参考能，不应该作为一个吸附物种写入最终形成能表。

---

## 步骤 E：建立工作能量 `E0`

```python
E0 = energy
```

这是非常值得记住的一行。

```text
energy
= 原始值，尽量不动

E0
= 工作变量，后面不断修改
```

可以把 `E0` 想成草稿纸。

以 OOH 为例：

```text
energy = -877.1433925 eV
E0     = -877.1433925 eV
```

---

## 步骤 F：如果是表面物种，先减裸板

```python
if site == '111':
    E0 -= references['111']
```

### 🐍 Python

```python
E0 -= x
```

等价于：

```python
E0 = E0 - x
```

### 🧪 OOH 实际计算

```text
E0
= -877.1433925 - (-863.25438)

= -13.8890125 eV
```

因此这一步不是最终形成能。

它只是完成：

> **完整“表面 + OOH”体系 → 去掉裸表面。**

### 状态快照

```text
开始：
表面 + OOH 总能
-877.1433925 eV

减去裸板：
- (-863.25438)

现在：
-13.8890125 eV

还需要：
继续扣 O、O、H 的参考能
```

---

## 步骤 G：把项目名称翻译成真正化学式

```python
formula = formula_map.get(
    name,
    name.replace('-', '')
)
```

### 🐍 Python

`dict.get(key, default)`：

```text
如果 key 存在 → 返回字典中的值
如果不存在   → 返回 default
```

### 🧪 为什么需要？

普通物种：

```text
OOH → OOH
H2O2 → H2O2
```

特殊物种：

```text
Ha   → H
Hb   → H
Hb-O → HO
```

这里 Python 并不知道 `Ha` 是“某个位点上的 H”。

因此必须人工告诉 ASE：

```python
formula_map['Ha'] = 'H'
```

---

## 步骤 H：解析元素组成

```python
try:
    composition = string2symbols(formula)
except Exception as exc:
    raise ValueError(...) from exc
```

### 🐍 Python

`try/except`：

```text
尝试执行 try
      │
      ├─ 成功 → 继续
      │
      └─ 出错 → 进入 except
```

`raise ... from exc` 会保留原始错误原因。

### 🧪 OOH 的结果

```text
formula = 'OOH'

string2symbols('OOH')
        ↓
composition = ['O', 'O', 'H']
```

所以 `composition` 不是一个“化学式字符串”了，而是一串具体原子。

---

## 步骤 I：逐原子扣参考能

```python
for atom in composition:
    if atom not in references:
        raise KeyError(...)
    E0 -= references[atom]
```

### 🐍 Python

♻ **复习：for 循环。**

这里会依次执行：

```text
atom = O
atom = O
atom = H
```

### 🧪 OOH 完整数值追踪

减裸板后：

```text
E0 = -13.8890125
```

第一次 O：

```text
E0 = -13.8890125 - (-4.37)
   = -9.5190125
```

第二次 O：

```text
E0 = -9.5190125 - (-4.37)
   = -5.1490125
```

H：

```text
E0 = -5.1490125 - (-1.11)
   = -4.0390125 eV
```

最终：

```python
formation_energies[key] = round(E0, 3)
```

得到：

```text
OOH_111 = -4.039 eV
```

### 一张图记住整个函数

```text
OOH_111 : -877.1434
        │
        ├── split
        │
        ├── name = OOH
        └── site = 111
                │
                ▼
       减裸板 -(-863.2544)
                │
                ▼
          -13.8890
                │
                ├── - Eref(O)
                ├── - Eref(O)
                └── - Eref(H)
                │
                ▼
           -4.0390
                │
                ▼
          round(...,3)
                │
                ▼
      formation_energies
```

---

# 1.3 真正调用函数

```python
formation_energies = get_formation_energies(
    abinitio_energies,
    ref_dict
)
```

### 🐍 Python

函数调用顺序：

```text
实参
  ↓
形参
  ↓
执行函数体
  ↓
return
  ↓
赋值给左侧变量
```

### 🧪 变量身份变化

```text
输入：
abinitio_energies
原始总能

      ↓ 函数

输出：
formation_energies
形成能
```

这两个字典绝不能混。

---

# 1.4 `frequency_dict`

```python
frequency_dict = {
    'H2O2_gas': [],
    ...
}
```

### 当前任务

> 🎯 给每个物种预留“振动频率”字段。

当前：

```text
[] = 暂时没有频率数据
```

所以它目前主要是为了满足 CatMAP 输入表结构，而不是在这里做振动热力学计算。

---

# 1.5 `make_input_file()`：把字典变成 CatMAP TXT

## 函数卡

| 项目 | 内容 |
|---|---|
| 输入 1 | `file_name`：输出文件名 |
| 输入 2 | `energy_dict`：这里传的是 `formation_energies` |
| 输入 3 | `frequencies`：这里传的是 `frequency_dict` |
| 输出 | 磁盘上的 TXT 文件 |
| 核心 Python | `join()`、列表、生成器、`with open()` |
| 核心项目逻辑 | 按 CatMAP TableParser 要求排列 6 列 |

### 参数对应关系

调用：

```python
make_input_file(
    file_name,
    formation_energies,
    frequency_dict
)
```

进入函数后：

```text
file_name   ← 'tini_energies.txt'
energy_dict ← formation_energies
frequencies ← frequency_dict
```

注意：

> 这里的 `energy_dict` 已经不再是原始总能，而是**形成能字典**。

这就是形参命名容易让人混淆的地方。

---

## 步骤 A：生成表头

```python
header = '\t'.join([
    'surface_name',
    'site_name',
    'species_name',
    'formation_energy',
    'frequencies',
    'reference',
])
```

### 🐍 Python

`\t` 是制表符。

`join()`：

```text
A + B + C
变成
A\tB\tC
```

### 🧪 输出含义

CatMAP 将收到：

```text
surface_name
site_name
species_name
formation_energy
frequencies
reference
```

这 6 列。

---

## 步骤 B：逐物种构造一行

```python
name, site = key.split('_', 1)

frequency = frequencies.get(key, [])
surface = None if site == 'gas' else 'tini'

outline = [
    surface,
    site,
    name,
    energy,
    frequency,
    'Input File Tutorial.',
]
```

### 🐍 Python

```python
A if condition else B
```

是条件表达式。

例如：

```text
site == gas
→ surface = None

site == 111
→ surface = tini
```

### 🧪 一行数据如何形成？

假设当前：

```text
key = OOH_111
energy = -4.039
```

则：

```text
name      = OOH
site      = 111
surface   = tini
frequency = []
```

最后形成：

```text
tini | 111 | OOH | -4.039 | [] | Input File Tutorial.
```

---

## 步骤 C：写文件

```python
with open(file_name, 'w', encoding='utf-8') as file:
    file.write(input_file)
```

### 🐍 Python

`with` 是上下文管理器。

优势：

```text
进入 with → 文件打开
离开 with → 自动关闭
```

比：

```python
f = open(...)
...
f.close()
```

更不容易忘记关闭文件。

### 状态快照

```text
formation_energies 仍在内存
        │
        ▼
make_input_file()
        │
        ▼
tini_energies.txt 已落盘
```

---

# 1.6 CatMAP 自检：不是正式跑微观动力学

```python
rxm = ReactionModel()
parser = TableParser(rxm)
parser.input_file = file_name
parser.parse()
```

## 当前任务

> 🎯 只验证刚才生成的 TXT **格式和物种定义是否能被 CatMAP 接受**。

这里：

```text
rxm
= 空/最小 ReactionModel

parser
= 负责读表的人
```

不是：

```text
正式完整微观动力学求解
```

### 对象关系

```text
ReactionModel()
      │
      ▼
     rxm
      │
      ▼
TableParser(rxm)
      │
      ▼
    parser
      │
      ├── input_file
      └── parse()
```

### 特殊物种

```python
rxm.species_definitions = {
    's': {'site_names': ['111']},
    'Ha_s': {'composition': {'H': 1}},
    'Hb_s': {'composition': {'H': 1}},
    'Hb-O_s': {'composition': {'H': 1, 'O': 1}},
}
```

这里的核心不是 Python，而是告诉 CatMAP：

```text
Ha_s   不是化学元素 Ha
Hb_s   不是化学元素 Hb
Hb-O_s = H1O1
```

---

# 1.7 本文件极简复现例子

只保留 OOH：

```python
energy = -877.1433925
slab = -863.25438

E0 = energy
E0 -= slab
E0 -= -4.37   # O
E0 -= -4.37   # O
E0 -= -1.11   # H

print(round(E0, 3))
```

输出应为：

```text
-4.039
```

如果这一小段完全理解，`get_formation_energies()` 的核心就已经理解。

---

# 2. `1-generating_input_file/generate_input1.py`

# 文件级固定 5 步

## 【1. 一句话总览】

这是上一个脚本的扩展版：仍然计算形成能，但额外输出 **TXT + CSV**，并强化 CatMAP 解析报错。

## 【2. 整体架构】

```text
原始能量
   ↓
形成能
   ↓
   ├── make_input_file() → TXT → CatMAP
   │
   └── make_csv_file()   → CSV → 人工检查/Excel
```

## 【3. 核心 Python 基础】

新的重点只有：

- `import csv`
- `csv.DictWriter`
- 关键字参数 `fieldnames=...`
- `isinstance()`
- 嵌套条件表达式
- `raise` 重新抛出异常

其余 `dict / for / split / try / with` 与上一文件相同。

## 【4. 关键变量 & 核心逻辑】

形成能计算逻辑不变。

新增逻辑：

```text
同一个 formation_energies
       │
       ├── TXT：给 CatMAP
       └── CSV：给人检查
```

## 【5. 记忆要点 + 易错坑】

> **TXT 是程序接口，CSV 是人工审查接口；两者必须来自同一套形成能。**

---

# 2.1 与上一文件重复的部分

以下部分已经学过，不再逐行重复：

```text
abinitio_energies
ref_dict
formula_map
get_formation_energies()
frequency_dict
make_input_file()
```

♻ **Python 复习：**

```text
字典 → 循环 → split → 扣参考 → return
```

真正值得新增学习的是 `make_csv_file()`。

---

# 2.2 `make_csv_file()`：CSV 输出器

## 函数卡

| 项目 | 内容 |
|---|---|
| 输入 | 文件名、形成能字典、频率字典 |
| 输出 | CSV 文件 |
| 核心 Python | `csv.DictWriter`、关键字参数、类型判断 |
| 目的 | 用 Excel/Origin 等方便地人工核查数据 |

### 创建 Writer

```python
with open(
    file_name,
    'w',
    newline='',
    encoding='utf-8'
) as csvfile:

    writer = csv.DictWriter(
        csvfile,
        fieldnames=fieldnames
    )

    writer.writeheader()
```

### 🐍 Python

`csv.DictWriter` 的核心思想：

> 不是按照第 1 列、第 2 列记位置，而是按“列名”写数据。

因此：

```python
writer.writerow({
    'surface_name': ...,
    'site_name': ...,
    'species_name': ...,
})
```

可读性比手工拼逗号更高。

---

## 频率格式化

```python
freq_str = ';'.join(
    str(int(freq))
    if isinstance(freq, (int, float)) and float(freq).is_integer()
    else str(freq)
    for freq in freq_list
) if freq_list else ''
```

这是本文件 Python 语法最密集的一段。

不要整句背，拆成 4 层：

```text
第 1 层：
for freq in freq_list
→ 一个一个取频率

第 2 层：
isinstance(freq, (int, float))
→ 是不是数值

第 3 层：
float(freq).is_integer()
→ 数值是不是整数值

第 4 层：
';'.join(...)
→ 多个频率用 ; 拼起来
```

例如：

```text
freq_list = [100.0, 250.5, 300]
```

得到：

```text
100;250.5;300
```

### 状态快照

```text
frequency_dict
    ↓
freq_list
    ↓ 格式化
freq_str
    ↓
CSV 一个单元格
```

---

# 2.3 CatMAP 解析失败为什么最后还要 `raise`

```python
except Exception as exc:
    print('CatMAP 解析测试失败：')
    print(exc)
    raise
```

### 🐍 Python

单独写：

```python
raise
```

表示重新抛出刚刚捕获的异常。

如果只：

```python
print(exc)
```

程序可能打印错误后仍然以“正常结束”退出。

加上：

```python
raise
```

就会：

```text
打印友好提示
+
保留 traceback
+
进程返回失败状态
```

这对调试和批处理都更可靠。

---

# 2.4 本文件状态快照

执行成功后：

```text
内存：
formation_energies

磁盘：
tini_energies.txt
tini_energies.csv

自检：
parser.parse() 成功
```

---

# 3. `2-creating_microkinetic_model/alkene_epoxidation.mkm`

`.mkm` 不是普通 `.py` 文件，但 CatMAP 按 Python 风格读取配置。因此它更适合当作：

> **“模型参数说明书”而不是“算法程序”。**

# 文件级固定 5 步

## 【1. 一句话总览】

告诉 CatMAP：**有哪些基元反应、有哪些表面、用什么描述符、什么反应条件、什么热力学近似、什么求解精度。**

## 【2. 整体架构】

```text
反应网络
rxn_expressions
      │
      ▼
表面 + 描述符
surface_names
descriptor_names
descriptor_ranges
      │
      ▼
反应条件
temperature
species_definitions
      │
      ▼
能量数据
input_file = energies.txt
      │
      ▼
热力学/标度
thermo_mode
scaling_constraint_dict
      │
      ▼
数值求解
precision / tolerance
```

## 【3. 核心 Python 基础】

这里主要复习：

- 变量赋值；
- 字符串；
- 列表；
- 嵌套列表；
- 字典；
- 嵌套字典；
- `None`；
- 科学计数法 `1e-20`。

没有复杂循环算法。

## 【4. 关键变量 & 核心逻辑】

这份文件真正需要记住的是 7 组配置：

```text
1. rxn_expressions       反应网络
2. surface_names         哪些催化剂
3. descriptor_names      横纵坐标是什么
4. descriptor_ranges     横纵坐标扫多大范围
5. species_definitions   边界条件/位点/组成
6. thermo/scaling        能量怎么处理
7. solver parameters     方程怎么求
```

## 【5. 记忆要点 + 易错坑】

> **机理决定“算什么”，描述符决定“在哪算”，条件决定“什么环境算”，solver 决定“怎么算出来”。**

---

# 3.1 `rxn_expressions`：反应网络

```python
rxn_expressions = [
    '* + H2O2_g <-> H2O2*',
    'H2O2* + * -> Ha* + OOH*',
    '* + C6H12_g <-> C6H12*',
    'OOH* + C6H12* + * <-> OOH-C6H12* + 2* -> C6H12O* + O* + Hb*',
    'C6H12O* <-> C6H12O_g + *',
    'Hb* + O* <-> Hb-O* + * -> OH* + *',
    'Ha* + OH* -> * + H2O*',
    'H2O* <-> H2O_g + *',
]
```

### 🐍 Python

Python 只把这些看成：

```text
一个 list
里面有 8 个 str
```

`<->`、`->`、`*` 的反应意义不是 Python 解释，而是 CatMAP 解析。

### 🧪 反应链

```text
H2O2 吸附
   ↓
H2O2 活化
   ↓
Ha* + OOH*
   ↓
烯烃吸附
   ↓
OOH* + C6H12*
   ↓
OOH-C6H12‡
   ↓
C6H12O* + O* + Hb*
   ↓
H 转移 / OH*
   ↓
H2O*
   ↓
H2O 脱附
```

### 如何读显式 TS

```text
初态 <-> 过渡态 -> 终态
```

例如：

```text
OOH* + C6H12* + *
        <->
OOH-C6H12* + 2*
        ->
C6H12O* + O* + Hb*
```

`OOH-C6H12*` 是显式过渡态名称。

---

# 3.2 `surface_names`：模型在哪些催化剂上建立

```python
surface_names = [
    'tife', 'timn', 'tihf', 'tire', 'timo', 'tiv', 'tizr',
    'tico', 'titi', 'tiw', 'tita', 'ticr', 'ti', 'tini'
]
```

### 参数卡

```text
变量：surface_names
类型：list[str]
含义：energies.txt 中的不同 Ti-M 催化剂/表面
用途：构建跨表面的能量-描述符关系
```

看到 `xxx_names`，优先理解成：

> “一组名字”，不是数值。

---

# 3.3 `descriptor_names` 与 `descriptor_ranges`

```python
descriptor_names = ['OOH_s', 'O_s']

descriptor_ranges = [
    [-7.0, -3.0],
    [-4.0,  2.5]
]

resolution = 20
```

## 当前任务

> 🎯 定义二维描述符地图。

配对关系必须按位置理解：

```text
descriptor_names[0] = OOH_s
descriptor_ranges[0] = [-7, -3]

descriptor_names[1] = O_s
descriptor_ranges[1] = [-4, 2.5]
```

所以：

```text
横/纵轴本质上是：
OOH_s
O_s
```

`resolution = 20` 表示每个方向使用 20 个网格点。

如果是二维网格，直觉上就是约：

```text
20 × 20
```

个描述符坐标点需要求解。

---

# 3.4 `temperature`

```python
temperature = 333.15
```

没有复杂 Python：

```text
变量名 = 浮点数
```

科研意义：

```text
333.15 K = 60 °C
```

它会进入速率常数和热力学相关计算。

---

# 3.5 `species_definitions`：边界条件数据库

这是 `.mkm` 最容易因参数多而迷失的一部分。

先把它分类，而不是逐行背。

```text
species_definitions
│
├── 储库/反应条件
│   ├── C6H12_g
│   ├── H2O2_g
│   ├── H2O_g
│   └── C6H12O_g
│
├── 位点
│   └── s
│
└── 特殊物种组成
    ├── Ha_s
    ├── Hb_s
    └── Hb-O_s
```

### 储库参数

```python
species_definitions['C6H12_g'] = {
    'pressure': 1.0
}

species_definitions['H2O2_g'] = {
    'pressure': 0.185
}

species_definitions['H2O_g'] = {
    'pressure': 0.815
}

species_definitions['C6H12O_g'] = {
    'pressure': 1e-20
}
```

当前模型把 `pressure` 作为液相体系的**有效活度/储库参数近似**，并不是严格气相分压。

记：

```text
C6H12   → 1.0
H2O2    → 0.185
H2O     → 0.815
产物    → 1e-20 ≈ 初始不存在
```

### 位点

```python
species_definitions['s'] = {
    'site_names': ['111'],
    'total': 1
}
```

```text
s
= CatMAP 内部抽象位点名

111
= 输入能量表中的实际 site_name
```

所以：

```text
s ↔ 111
```

### 非标准物种

```python
species_definitions['Ha_s'] = {
    'composition': {'H': 1}
}
```

仍然是在解决同一个问题：

> CatMAP/ASE 不知道 Ha 是 H，需要人工指定元素组成。

---

# 3.6 数据文件

```python
data_file = 'alkene_epoxidation.pkl'
input_file = 'energies.txt'
```

不要混：

```text
input_file
= 输入
= energies.txt

data_file
= CatMAP 运行数据/缓存
= alkene_epoxidation.pkl
```

---

# 3.7 热力学模式

```python
gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
pressure_mode = 'static'
```

快速记忆：

```text
frozen_gas
→ 当前储库物种不额外采用完整气相热修正

frozen_adsorbate
→ 吸附物采用冻结近似

static
→ 边界 pressure/activity 不随计算动态改变
```

---

# 3.8 `scaling_constraint_dict`

```python
scaling_constraint_dict = {
    'OOH_s': ['+', 0, None],
    'O_s':   [0, '+', None],
}
```

## 当前任务

> 🎯 约束两个被选作 descriptor 的物种与两个描述符方向对应。

不要把这一行当普通 Python 算法。

Python 只负责保存：

```text
key → list
```

真正的 `+ / 0 / None` 由 CatMAP scaler 解释。

直觉上：

```text
OOH_s
主要沿 descriptor 1

O_s
主要沿 descriptor 2
```

---

# 3.9 Solver 参数

```python
decimal_precision = 100
tolerance = 1e-50
max_rootfinding_iterations = 100
max_bisections = 3
```

把它们按功能分组：

```text
decimal_precision
→ 数字算多精细

tolerance
→ 多小才算收敛

max_rootfinding_iterations
→ 最多迭代多少次

max_bisections
→ 最多辅助二分多少次
```

### 状态快照

`.mkm` 执行到末尾时，还没有画图。

它只是完成：

```text
✓ 反应网络
✓ 催化剂集合
✓ 描述符
✓ 条件
✓ 能量来源
✓ scaler
✓ solver
```

真正求解发生在：

```python
model.run()
```

---

# 4. `2-creating_microkinetic_model/test.ipynb`

# 文件级固定 5 步

## 【1. 一句话总览】

读取 `.mkm` 或 `.log`，真正运行 CatMAP，并把二维微观动力学结果绘图、导出和检查。

## 【2. 整体架构】

```text
Cell 1
.mkm → ReactionModel → run → VectorMap / ScalingAnalysis

Cell 2
.log → ReactionModel → 结果 map → table.txt

Cell 3
.log → coverage_map → 打印覆盖度
```

## 【3. 核心 Python 基础】

新重点：

- 类实例化；
- 对象属性；
- 对象方法；
- `+=` 扩展列表；
- `glob()` 搜文件；
- `len()`；
- 列表索引；
- 自定义辅助函数；
- 嵌套循环；
- 列表推导式；
- `getattr()` 动态属性。

## 【4. 关键变量 & 核心逻辑】

```text
model
= 数据和计算结果的中心对象

vm
= 从 model 拿数据并画二维图

sa
= 从 model 做 scaling analysis
```

## 【5. 记忆要点 + 易错坑】

> **model 负责算，vm 负责画，log 负责恢复。**

当前 Notebook 有几处旧示例残留，需要特别注意：

- `H reactivity` / `C2H5 reactivity` 与当前 `OOH_s/O_s` 不一致；
- `C2H5_s` 当前网络不存在；
- `MgO` 变量名与 Ti-M 项目不一致；
- Cell 2 若独立运行，`model.output_variables += ...` 的位置存在问题；
- 未处理没有 `.log` 文件的情况；
- `output_labels['production_rate']` 被写死。

---

# 4.1 Cell 1：正式运行模型

## 当前任务

> 🎯 从 `.mkm` 构建模型，然后调用 `run()` 求稳态。

```python
from catmap import ReactionModel

mkm_file = 'alkene_epoxidation.mkm'
model = ReactionModel(setup_file=mkm_file)

model.output_variables += [
    'production_rate',
    'rate',
    'rate_control',
    'coverage',
    'selectivity_control',
    'rxn_order'
]

model.run()
```

## 参数卡

```text
mkm_file
= 配置文件路径

model
= 完整 ReactionModel 对象

output_variables
= 希望 CatMAP 额外保存哪些结果

run()
= 真正开始数值求解
```

### 🐍 `+=` 在列表中的意义

```python
model.output_variables += ['coverage']
```

可以理解为：

```text
原来的 output_variables
+
新的输出名字
=
扩展后的列表
```

### 🧪 状态变化

运行前：

```text
model
有机理/参数
没有完整网格结果
```

执行：

```python
model.run()
```

之后：

```text
model
├── rate_map
├── production_rate_map
├── coverage_map
├── rate_control_map
└── ...
```

因此 `_map` 可以记成：

> **“某个输出变量在描述符网格上的结果集合”。**

---

# 4.2 `VectorMap`：不是重新计算，是画已有结果

```python
from catmap import analyze

vm = analyze.VectorMap(model)
vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-25
vm.max = 1e3
vm.plot(save='rate.pdf')
```

## 对象关系

```text
model
  │
  │ 已经有结果
  ▼
VectorMap(model)
  │
  ▼
 vm
  │
  ├── plot_variable
  ├── log_scale
  ├── min/max
  └── plot()
```

看到：

```python
vm.plot_variable = 'coverage'
```

不要理解成“计算 coverage”。

它只是：

> 告诉绘图器下一张图画 `coverage`。

---

# 4.3 同一个 `vm` 为什么不断改属性？

例如：

```python
vm.plot_variable = 'rate'
vm.plot(...)

vm.plot_variable = 'production_rate'
vm.plot(...)

vm.plot_variable = 'coverage'
vm.plot(...)
```

这是典型的**对象状态复用**。

```text
同一个 vm
  │
  ├── 改成 rate → 画
  ├── 改成 production_rate → 画
  └── 改成 coverage → 画
```

所以阅读对象代码时，必须问：

> **这个属性是新建了一个对象，还是只修改现有对象的状态？**

这里属于后者。

---

# 4.4 当前两个旧残留

```python
vm.descriptor_labels = [
    'H reactivity [eV]',
    'C2H5 reactivity [eV]'
]
```

但当前 `.mkm`：

```python
descriptor_names = ['OOH_s', 'O_s']
```

因此标签不一致。

另外：

```python
vm.include_labels = ['C2H5_s']
```

当前网络也没有 `C2H5_s`。

理解代码时，应区分：

```text
Python 语法正确
≠
项目逻辑正确
```

---

# 4.5 ScalingAnalysis

```python
sa = analyze.ScalingAnalysis(model)
sa.plot(save='scaling.pdf')
```

变量记忆：

```text
sa
= Scaling Analysis
```

它不负责主求解，而是对 `model` 中的能量/标度关系做分析和绘图。

---

# 4.6 Cell 2：从 `.log` 恢复已有模型

```python
from glob import glob

logfile = glob('*.log')
```

### 🐍 Python

`glob('*.log')` 返回一个列表：

```text
[
  'alkene_epoxidation.log',
  ...
]
```

`*` 表示任意文件名。

因此：

```python
len(logfile)
```

表示找到多少个日志。

---

## 这里应当有三个分支

当前代码主要检查：

```python
if len(logfile) > 1:
```

完整逻辑最好理解成：

```text
len == 0
→ 没找到 log

len == 1
→ 正常

len > 1
→ 不知道该选哪个
```

这是典型的边界情况思维。

---

# 4.7 `flatten_2d()`：二维结果摊平

```python
def flatten_2d(output):
    flat = []
    for x in output:
        flat += x
    return flat
```

## 输入输出

例如：

```text
输入：
[
  [1, 2],
  [3, 4]
]

输出：
[1, 2, 3, 4]
```

### 🐍 Python

列表中：

```python
flat += x
```

相当于把 `x` 中的元素一个个追加进去。

它和：

```python
flat.append(x)
```

不同：

```text
+=      → [1,2,3,4]
append  → [[1,2],[3,4]]
```

---

# 4.8 `getattr()`：根据字符串动态选结果

```python
for pt, output in getattr(
    model,
    output_variable + '_map'
):
```

这是 Cell 2 最值得学习的 Python 技巧。

假设：

```python
output_variable = 'production_rate'
```

那么：

```python
output_variable + '_map'
```

得到：

```text
'production_rate_map'
```

然后：

```python
getattr(model, 'production_rate_map')
```

等价于：

```python
model.production_rate_map
```

因此它实现：

```text
变量名字决定读取哪个属性
```

避免写很多：

```python
if output_variable == 'rate':
    ...
elif output_variable == 'coverage':
    ...
```

---

# 4.9 `pt` 与 `output` 是什么？

```python
for pt, output in model.production_rate_map:
```

每次循环可以理解成：

```text
pt
= 当前描述符坐标
例如 [E_OOH, E_O]

output
= 这个坐标下的一整组 production rates
```

所以数据结构是：

```text
descriptor point
       │
       ▼
microkinetic output
```

最终写表：

```text
descriptor-OOH_s
descriptor-O_s
各个 production rate
```

---

# 4.10 Cell 3：覆盖度检查

```python
labels = model.output_labels['coverage']

for MgO, cvg in model.coverage_map:
    print('descriptors', MgO)
    print('intermediates', labels)
    print('coverages', [float(c) for c in cvg])
```

建议读成：

```text
labels
= 每个覆盖度数值对应哪个物种

MgO
= 实际上是描述符坐标
  只是变量名沿用了旧项目

cvg
= coverage vector
```

所以更合适的理解名称其实是：

```python
for descriptor_point, coverage_vector in model.coverage_map:
```

### Python 与项目语义的区别

Python 完全不在乎变量叫：

```text
MgO
abc
x
descriptor_point
```

但人会在乎。

所以科研代码应优先使用有物理意义的变量名。

---

# 5. 贯穿项目的“参数不遗忘”阅读法

以后看到陌生变量，不要试图硬记。固定问四个问题：

```text
1. 它是什么类型？
2. 它从哪里来？
3. 此刻里面装的是什么？
4. 下一步谁会使用它？
```

例如：

```python
composition = string2symbols(formula)
```

立刻做四问：

```text
类型？
→ list

从哪里来？
→ ASE 函数返回

里面是什么？
→ ['O','O','H']

下一步谁用？
→ for atom in composition
```

这比记“composition 中文是组成”有效得多。

---

# 6. Python 基础知识复习路线

不要按照 Python 教材顺序重新学一遍。

直接按照当前项目出现频率：

```text
第 1 级
变量
dict
list
for
if

第 2 级
函数
形参/实参
return
split/get/join

第 3 级
异常
with open
模块 import

第 4 级
类与对象
对象属性
对象方法

第 5 级
列表推导式
生成器表达式
getattr
嵌套数据结构
```

学完第 1~4 级，已经可以看懂本项目绝大多数代码。

---

# 7. 一页式总复习

```text
【能量层】

abinitio_energies
原始总能
      │
      ▼
get_formation_energies()
      │
      ▼
formation_energies
形成能


【文件层】

formation_energies
      │
      ├── make_input_file()
      │       ↓
      │      TXT
      │
      └── make_csv_file()
              ↓
             CSV


【模型层】

energies.txt
     +
alkene_epoxidation.mkm
      │
      ▼
ReactionModel
      │
      ▼
model.run()


【结果层】

model
├── rate_map
├── production_rate_map
├── coverage_map
└── rate_control_map
      │
      ▼
VectorMap / table / print
```

---

# 8. 当前代码中最值得优先修改/检查的地方

1. `surface = None if site == 'gas' else 'tini'` 将表面名硬编码成 `tini`。若后续自动批量处理 TiFe、TiMn、TiHf 等，应把表面名变成函数参数。
2. `test.ipynb` 的 `H reactivity`、`C2H5 reactivity` 与当前 `OOH_s / O_s` 描述符不一致。
3. `vm.include_labels = ['C2H5_s']` 与当前反应网络不一致。
4. `import page`、变量名 `MgO` 属于旧示例痕迹。
5. Notebook Cell 2 应先创建/恢复 `model`，再访问 `model.output_variables`。
6. `glob('*.log')` 应同时处理 0 个、1 个、多个日志文件。
7. `model.output_labels['production_rate']` 若目标是通用导出，应考虑使用 `model.output_labels[output_variable]`。

---

# 9. 最后只记 6 句话

```text
1. dict 是这个项目最核心的数据容器。
2. energy 是原始值，E0 是计算中的工作值。
3. name 是项目名字，formula 是化学式，composition 是元素列表。
4. formation_energies 是 DFT 总能经过参考态换算后的结果。
5. .mkm 主要负责“配置”，model.run() 才是真正求解。
6. model 负责算，vm 负责画，*_map 保存网格结果。
```

以后新增代码文件时，继续沿用：

> **文件级固定 5 步 + 单文件 D 双轨科研伴读法。**
