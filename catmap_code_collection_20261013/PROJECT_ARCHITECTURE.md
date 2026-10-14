# CatMAP 项目归档架构约定

统一归档根目录：`H:\catmap\catmap\catmap_code_collection_20261013`

后续 CatMAP 相关模型、脚本、计算结果和分析图优先放入该归档根目录，并按工作性质分区；禁止把新结果散放在根目录，禁止覆盖原始数据、已验收结果或历史快照。

## 目录约定

```text
catmap_code_collection_20261013/
├── catmap_official_run/         # 官方 CatMAP 输入、模型文件和官方运行结果
├── models/                      # 经阶段标识的模型配置/历史模型
├── scripts/                     # 项目级审计、转换和通用脚本
├── tests/                       # 回归测试
├── v2_global_reference/         # V2 全局参考审计
├── analyses/                    # 派生分析；每个研究问题独立子目录
│   └── <analysis_name>/
│       ├── data/                # 输入快照、匹配表和来源审计
│       ├── figures/             # PNG/PDF 图片
│       ├── origin_exports/      # Origin/其他作图软件可读的网格与点表
│       ├── scripts/             # 该分析专属的可复现脚本
│       └── reports/             # 方法说明、结果说明和限制
└── PROJECT_ARCHITECTURE.md      # 本架构约定
```

## 数据与结果规则

- 原始 Excel、DFT/NEB 数据和正式基线保持只读；派生数据必须保留来源说明及必要的哈希。
- 官方 CatMAP 输入/运行与经验统计分析分开放置，经验响应面不得标成 CatMAP 连续计算结果。
- 目录名应标明研究阶段/用途；新一轮输出应新建目录或使用新文件名，不覆盖历史结果。
- 将 C 盘工作区内容归档到此处时，默认复制而不移动；完成后校验文件数、大小及 SHA256。
- 运行脚本若仍引用 C 盘绝对路径，应在报告中说明；没有完成路径迁移和复现验证前，不宣称归档副本可独立运行。

## 本次新增分析

O*/Ha*/Hb* 历史形成能与严格 CatMAP TOF 的经验二维图、全矩形 Origin 网格及生成脚本归档至：

`analyses/empirical_O_H_descriptor_analysis/`

13个位点 CatMAP 微观动力学机制分析（稳态覆盖度、六步通量、有效势垒、官方DRC及TS1–TS2网格覆盖度复核）归档至：

`analyses/microkinetic_mechanism_13sites/`

后续 CatMAP 相关代码和结果继续遵循本文件的功能分区：官方模型/运行放入 `catmap_official_run/`，通用项目脚本与测试分别放入 `scripts/`、`tests/`，研究问题专属的派生数据、脚本、图和报告放入对应的 `analyses/<analysis_name>/` 子目录。归档默认复制，不覆盖原始输入、正式基线或历史分析。
