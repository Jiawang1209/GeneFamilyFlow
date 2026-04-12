# GeneFamilyFlow

多物种基因家族自动化分析流程。Snakemake 编排，Python 处理数据，R 做统计和可视化。

## 技术栈

| 层级 | 工具 |
|------|------|
| 流程管理 | Snakemake |
| 数据处理 | Python 3.11+ |
| 统计可视化 | R 4.x |
| 序列分析 | HMMER、BLAST、seqkit、MUSCLE、IQ-TREE、MEME |
| 共线性 | JCVI、MCScanX |

## 目录结构

```
GeneFamilyFlow/
├── Snakefile                  # 主流程（待创建）
├── config/default_config.yaml # 全流程参数
├── envs/genefamily.yaml       # Conda 环境
├── scripts/                   # Python 工具（含单元测试）
├── R/                         # R 可视化脚本（optparse 接口）
├── tests/                     # pytest 测试
├── example/                   # 真实案例数据和手工分析记录
├── docs/                      # 详细文档（setup/faq/contributing）
├── output/                    # 流程输出（运行后生成）
└── work/                      # 中间文件（运行后生成）
```

## 代码规范

### Python（`scripts/`）

- 完整类型注解（`str | Path`、`list[T]`、`Iterator[T]`）
- 数据类用 `@dataclass(frozen=True)`；不修改对象，返回新对象
- 所有脚本提供 argparse CLI，入口 `main(argv) -> int`
- 自定义异常继承清晰的基类（如 `FastaFormatError(ValueError)`）
- 大文件用生成器流式处理
- **每个脚本必须有对应 `tests/test_*.py`**

参考已有实现：`scripts/filter_longest_transcript.py`、`scripts/read_protein_fasta.py`。

### R（`R/`）

- 所有脚本用 `optparse` 接收参数，便于 Snakemake 集成
- tidyverse 风格 + ggplot2
- 输出 PDF（矢量）到 `--outdir`
- 必需参数缺失时 `print_help()` + `quit(status=1)`

参考已有实现：`R/05_genefamily_info.R`、`R/06_tree.R`。

### Snakemake

- 每个规则对应一个分析步骤，命名 `stepN_description`
- 输入/输出路径从 `config` 读取，不硬编码
- 每个规则声明 `conda: "envs/genefamily.yaml"`
- 大文件用 `temp()`，日志用 `log:`
- 编写后先 dry-run：`snakemake -n -p`

## 常用命令

```bash
# 环境
conda env create -f envs/genefamily.yaml -n genefamily && conda activate genefamily

# 测试
pytest tests/ -v

# 流程
snakemake --configfile config/default_config.yaml -n -p           # dry-run
snakemake --configfile config/default_config.yaml -j 8 --use-conda # 运行
snakemake --configfile config/default_config.yaml --until step5_genefamily_info -j 4  # 运行到某步
```

## 添加新分析步骤

1. `scripts/` 加 Python 工具（带类型注解和 CLI）
2. `tests/` 加单元测试，跑通 `pytest`
3. `R/` 加可视化脚本（optparse 接口）
4. `Snakefile` 加规则，更新 `rule all` 输入
5. `config/default_config.yaml` 加该步骤参数
6. `snakemake -n` 验证 DAG，再实际运行

详细教程见 `docs/contributing.md`。

## 详细文档

- `docs/setup.md` — 环境搭建与故障排除
- `docs/faq.md` — 常见问题
- `docs/contributing.md` — 新增步骤完整示例
- `example/GeneFamilyFlow.md` — 11 步分析工作流说明

## 协作

- 规划/审查：使用 `/ask codex` 进行计划和代码审查（见 `~/.claude/CLAUDE.md` 的 Peer Review Framework）
- 多 Agent 规则：`.agentbridge/collaboration.md`
