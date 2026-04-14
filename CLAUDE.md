# GeneFamilyFlow

多物种基因家族自动化分析流程。Snakemake 编排，Python 处理数据，R 做统计和可视化。

## 技术栈

| 层级 | 工具 |
|------|------|
| 流程管理 | Snakemake 8+ |
| 数据处理 | Python 3.11+ |
| 统计可视化 | R 4.3+ |
| 序列分析 | HMMER、BLAST、seqkit、MUSCLE、ClustalW、pfam_scan.pl |
| 建树 | IQ-TREE 或 FastTree（`step6.tree_tool` 切换） |
| Motif | MEME |
| 共线性 | JCVI、MCScanX（支持 `precomputed` 模式） |
| 启动子 | bedtools + PlantCARE（在线） |

## 目录结构

```
GeneFamilyFlow/
├── Snakefile                  # 主流程（36-55 jobs）
├── config/default_config.yaml # 全流程参数
├── envs/                      # Conda 环境定义
├── scripts/                   # Python 工具（含单元测试）
├── R/                         # R 可视化脚本（optparse CLI）
├── tests/                     # pytest 测试（361 tests, 94% cov）
├── example/                   # 案例数据
├── output/                    # 流程输出（运行后生成）
└── work/                      # 中间文件（运行后生成）
```

## 流程架构

- **Steps 1-4**：数据准备 → 多 domain HMM 搜索 → BLAST → 合并鉴定
  - `{species}` × `{domain}` 双 wildcard 交叉
  - per-domain 独立搜索 → union → pfam_scan 验证全部 domain
- **Step 5**：基因家族理化性质统计
- **Step 6**：建树（IQ-TREE 或 FastTree）
- **Step 7**：Motif + 基因结构复合图
- **Steps 8-9**：共线性（JCVI/MCScanX，`precomputed` 可选）
- **Steps 10-11**：启动子分析 + PPI 网络
- **Steps 13-14**：可选（GO/KEGG、qRT-PCR）

## 代码规范

### Python（`scripts/`）

- 完整类型注解（`str | Path`、`list[T]`、`Iterator[T]`）
- 数据类用 `@dataclass(frozen=True)`；不修改对象，返回新对象
- 所有脚本提供 argparse CLI，入口 `main(argv) -> int`
- 大文件用生成器流式处理
- **每个脚本必须有对应 `tests/test_*.py`**

### R（`R/`）

- 所有脚本用 `optparse` 接收参数，便于 Snakemake 集成
- tidyverse 风格 + ggplot2，输出 PDF 到 `--outdir`
- 必需参数缺失时 `stop()` 报错

### Snakemake

- 每个规则对应一个分析步骤，命名 `stepN_description`
- 输入/输出路径从 `config` 读取，不硬编码
- 条件规则用 `if/else`（如 tree_tool、precomputed）
- 大文件用 `temp()`，日志用 `log:`
- 编写后先 dry-run：`snakemake -n -p`

## 常用命令

```bash
# 环境
conda env create -f envs/genefamily.yaml
conda activate GeneFamilyFlow

# 测试
pytest tests/ -v

# 流程
snakemake --configfile config/default_config.yaml -n -p           # dry-run
snakemake --configfile config/default_config.yaml -j 8 --use-conda # 运行
snakemake --configfile config/default_config.yaml --until step06_tree_build -j 4
```

## 添加新分析步骤

1. `scripts/` 加 Python 工具（带类型注解和 CLI）
2. `tests/` 加单元测试，跑通 `pytest`
3. `R/` 加可视化脚本（optparse 接口）
4. `Snakefile` 加规则，更新 `rule all` 输入
5. `config/default_config.yaml` 加该步骤参数
6. `snakemake -n` 验证 DAG，再实际运行
