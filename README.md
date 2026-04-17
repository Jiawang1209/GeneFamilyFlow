# GeneFamilyFlow

[![tests](https://github.com/Jiawang1209/GeneFamilyFlow/actions/workflows/test.yml/badge.svg)](https://github.com/Jiawang1209/GeneFamilyFlow/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://www.python.org/)

多物种基因家族自动化分析流程。Snakemake 编排 + Python 数据处理 + R 统计可视化，从蛋白/基因组 FASTA 一路跑到建树、共线性、启动子、PPI 出图。

整合了 HMM 搜索、BLAST、Pfam 验证、建树、MEME motif、JCVI/MCScanX 共线性、FIMO 启动子扫描、PPI 可视化，全部参数集中在一个 YAML 里。

---

## 目录

- [功能亮点](#功能亮点)
- [流程总览](#流程总览)
- [快速开始（5 分钟跑完案例）](#快速开始5-分钟跑完案例)
- [新机器从零搭环境](#新机器从零搭环境)
- [用自己的数据跑流程](#用自己的数据跑流程)
  - [1. 数据放哪里](#1-数据放哪里)
  - [2. 写配置文件](#2-写配置文件)
  - [3. 一键运行](#3-一键运行)
  - [4. 结果在哪里找](#4-结果在哪里找)
- [分步/按需运行](#分步按需运行)
- [ARM Mac 特别说明（Step 10 Docker）](#arm-mac-特别说明step-10-docker)
- [测试与开发](#测试与开发)
- [文档索引](#文档索引)
- [引用](#引用)

---

## 功能亮点

- **多 Pfam 结构域支持** — 每个结构域独立 HMM/BLAST，Step 4 按 `any`/`all` 合并
- **多物种并行** — 通过 `{species}` wildcard 自动 fan-out
- **建树二选一** — IQ-TREE（准）或 FastTree（快），配置切换
- **自动亚家族聚类** — 自顶向下单系拆分，圆形 + 矩形两种布局同时出图
- **Precomputed 模式** — JCVI/MCScanX 可跳过重算，直接出图
- **Step 10 离线默认** — 内置 JASPAR 2024 plants 束（805 个 PFM），FIMO 扫描，无需 PlantCARE 上传
- **单一配置文件** — 所有参数在 `config/default_config.yaml`
- **测试覆盖** — 389 个 pytest 用例，97% 覆盖率

## 流程总览

| Step | 作用 | 主要工具 |
|------|------|---------|
| 1 | 数据准备（最长转录本过滤） | seqkit |
| 2 | HMM 搜索（每个结构域 2 轮） | HMMER + ClustalW |
| 3 | BLAST 搜索（每个结构域） | BLAST+ |
| 4 | 基因家族鉴定（合并 + Pfam 验证） | pfam_scan.pl |
| 5 | 基因家族理化性质 | R (Biostrings, Peptides) |
| 6 | 建树 + 亚家族聚类 | MUSCLE + IQ-TREE/FastTree + ggtree |
| 7 | Motif + 基因结构复合图 | MEME + R (ggtree, gggenes) |
| 8 | JCVI 共线性 + Ka/Ks | JCVI + R |
| 9 | MCScanX 共线性 + Circos | MCScanX + R (circlize) |
| 10 | 启动子顺式元件 | bedtools + FIMO/JASPAR + R |
| 11 | PPI 网络 | R (ggNetView) |
| 13 *(可选)* | GO/KEGG 富集 | eggNOG-mapper + clusterProfiler |
| 14 *(可选)* | qRT-PCR 柱状图 | R (ggpubr) |

详细架构与 I/O 见 [CLAUDE.md](CLAUDE.md) 和 [docs/tutorial.md](docs/tutorial.md)。

---

## 快速开始（5 分钟跑完案例）

仓库自带 NPF 基因家族（*Sorghum bicolor* / *Arabidopsis thaliana* / *Oryza sativa*）的完整样例数据，Step 8/9 默认 precomputed，可直接跑通。

```bash
# 1. clone 仓库
git clone https://github.com/Jiawang1209/GeneFamilyFlow.git
cd GeneFamilyFlow

# 2. 建 conda 环境（~5 min）
conda env create -f envs/genefamily.yaml
conda activate GeneFamilyFlow

# 3. 先 dry-run 看看 DAG
snakemake --configfile config/default_config.yaml -n -p

# 4. 跑到 Step 5（不需要 Pfam 数据库、不需要 MEME）
snakemake --configfile config/default_config.yaml --until step05_genefamily_info -j 4 --use-conda

# 5. 看结果
open output/05_genefamily_info/Gene_Information.xlsx
```

跑全流程（含 MEME/JCVI/MCScanX/FIMO/PPI）需要先装 [Pfam 数据库](#步骤-4-可选-下载-pfam-数据库仅-step-4-需要) 和几个外部工具，详见下文。

---

## 新机器从零搭环境

> 目标系统：Linux（推荐）或 macOS（Intel / Apple Silicon 都行，M 系列芯片 Step 10 要用 Docker，见下文）

### 步骤 1. 装 Conda/Mamba

```bash
# 推荐 miniforge（自带 conda-forge）
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
bash Miniforge3-*.sh
```

### 步骤 2. 克隆仓库

```bash
git clone https://github.com/Jiawang1209/GeneFamilyFlow.git
cd GeneFamilyFlow
```

### 步骤 3. 建 Conda 环境

```bash
# 完整环境（所有 Python + R + 生信工具）
conda env create -f envs/genefamily.yaml
conda activate GeneFamilyFlow
```

完整环境包含：snakemake、HMMER、BLAST+、seqkit、MUSCLE、ClustalW、IQ-TREE、FastTree、MEME、FIMO、bedtools、pfam_scan.pl、R + tidyverse/ggtree/gggenes/circlize 等。

如果只想做开发/测试，不需要 R 生信工具：

```bash
conda env create -f envs/genefamily-minimal.yaml
conda activate GeneFamilyFlow
```

### 步骤 4. （可选）下载 Pfam 数据库（仅 Step 4 需要）

Step 4 用 `pfam_scan.pl` 验证结构域，需要本地 Pfam-A 库：

```bash
mkdir -p ~/pfam_db && cd ~/pfam_db
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip *.gz
hmmpress Pfam-A.hmm
```

然后把路径写进 config：

```yaml
step4:
  pfam_database_dir: "/Users/YOURNAME/pfam_db"
```

如果暂时不想装 Pfam 库，把 `step4.use_pfam_scan` 设为 `false` 即可跳过。

### 步骤 5. （可选）装 Conda 没有的外部工具

| 工具 | 用途 | 安装 |
|------|------|------|
| MCScanX | Step 9 共线性 | [GitHub](https://github.com/wyp1125/MCScanX) 源码编译 |
| KaKs_Calculator | Step 8/9 Ka/Ks | [SourceForge](https://sourceforge.net/projects/kakscalculator2/) |
| ggNetView | Step 11 PPI | `Rscript -e 'install.packages("ggNetView")'` |

Step 8/9 可以用 `precomputed: true` 跳过重算，直接拿现成结果出图 → 这两个工具可以暂时不装。

### 步骤 6. 验证环境

```bash
# 测试全绿 = Python 侧 OK
pytest tests/ -v

# dry-run 无错 = DAG 能解析、配置合法
snakemake --configfile config/default_config.yaml -n -p
```

---

## 用自己的数据跑流程

### 1. 数据放哪里

默认数据目录是 `example/1.database/`，但你可以新建一个自己的目录（推荐）。命名规范是 `{species}.{type}.fasta/gff3`：

```
your_data/
├── 1.database/                     # database_dir
│   ├── Species1.pep.fasta          # [必需] 蛋白序列
│   ├── Species1.gff3               # [step 5/8/9/10 需要] 基因注释
│   ├── Species1.cds.fasta          # [step 8 non-precomputed 需要] CDS
│   ├── Species1.genome.fasta       # [step 10 compute_promoter 需要] 基因组
│   ├── Species2.pep.fasta
│   └── ...
├── 2.hmmsearch/
│   └── PFxxxxx.hmm                 # [必需] 每个 Pfam 结构域的 HMM 模型
├── 3.blast/references/
│   ├── Athaliana.domains.tsv       # [必需] 2 列 TSV: gene_id<TAB>pfam_id
│   ├── Athaliana.pep.fasta         # [必需] 参考物种全蛋白组
│   └── ...
└── 10.promoter/
    └── species_10.bed              # [step 5/9/10 需要] 基因坐标 BED
                                    # 6 列: Chr Start End ID Info Strand
```

**最少需要的文件**（只跑到 Step 5）：

- 每个物种的 `{species}.pep.fasta`
- 每个 Pfam 结构域的 `.hmm` 文件
- 至少一个参考物种的 `.domains.tsv` + `.pep.fasta`（用来构建 BLAST seed）
- `species_10.bed` 基因坐标表

**注意：**
- FASTA header 第一段（到第一个空格前）必须是 gene ID
- `domains.tsv` 是两列 TSV：`gene_id<TAB>PF00854`，一行一个命中
- 可以把数据放在仓库外任何地方，config 里写绝对路径即可

### 2. 写配置文件

复制默认配置改一份：

```bash
cp config/default_config.yaml config/my_config.yaml
```

**最少要改的字段：**

```yaml
# 要分析的物种列表
species:
  - Species1
  - Species2

# 主物种（Step 6-10 以它为核心）
target_species: "Species1"

# 物种 → 基因 ID 前缀的映射（Step 5/6/8 用来按物种分颜色）
species_id_prefix:
  Species1: "Sp1."
  Species2: "AT"

# 数据路径
database_dir: "/path/to/your_data/1.database"
output_dir: "output"
work_dir: "work"

# Step 2：每个 Pfam 结构域一条
step2:
  domains:
    - pfam_id: "PF00854"
      pfam_hmm_file: "/path/to/your_data/2.hmmsearch/PF00854.hmm"
      seed_references:
        - Athaliana                 # 必须在 step3.reference_species 里

# Step 3：参考物种（做 BLAST seed 用）
step3:
  reference_species:
    Athaliana:
      domains_table: "/path/to/your_data/3.blast/references/Athaliana.domains.tsv"
      proteome:      "/path/to/your_data/3.blast/references/Athaliana.pep.fasta"

# Step 4：Pfam 库路径（或关掉 pfam_scan）
step4:
  use_pfam_scan: true
  pfam_database_dir: "/Users/YOURNAME/pfam_db"

# Step 5：基因坐标 BED
step5:
  gene_bed_file: "/path/to/your_data/10.promoter/species_10.bed"

# Step 6：建树工具 + 亚家族数
step6:
  tree_tool: "iqtree"               # 或 "fasttree"
  auto_n_subfamilies: 8             # 自动切分 8 个单系亚家族
  tree_layout: "both"               # "circular" / "rectangular" / "both"

# Step 8/9：有现成结果就用 precomputed
step8:
  precomputed: true                 # 跳过 JCVI 重算
step9:
  precomputed: true
  species_config:                   # 每个物种的染色体数 + genome .fai
    Species1:
      n_chromosomes: 10
      genome_length_file: "..."
      npf_id_file: "..."

# Step 10：启动子扫描方法
step10:
  compute_promoter: false           # true = 从基因组提启动子；false = 用现成 fasta
  scan_method: "jaspar"             # 推荐：离线 FIMO + JASPAR plants
  upstream_distance: 2000
```

其余可选步骤（13 GO/KEGG、14 qRT-PCR）默认 `enabled: false`，要用再开。

完整参数说明见 [docs/configuration.md](docs/configuration.md)。

### 3. 一键运行

```bash
# 先 dry-run 看 DAG 有没有问题
snakemake --configfile config/my_config.yaml -n -p

# 正式跑（8 并行）
snakemake --configfile config/my_config.yaml -j 8 --use-conda

# 失败了想从断点继续？Snakemake 默认就支持，直接再跑一次即可
# 想强制重跑某一步：
snakemake --configfile config/my_config.yaml --forcerun step06_tree_build -j 4 --use-conda
```

### 4. 结果在哪里找

所有最终图表都在 `output/`，中间文件在 `work/`，日志在 `logs/`。

```
output/
├── 04_identification/
│   └── identify.ID.clean.fa                # 最终基因家族序列（所有物种）
├── 05_genefamily_info/
│   ├── Gene_Information.xlsx               # 理化性质表（长度/MW/pI/亚细胞定位等）
│   └── Gene_Information.csv
├── 06_tree/
│   ├── phylogenetic_tree_circular.pdf      # 圆形建树 + 亚家族高亮
│   └── phylogenetic_tree_rectangular.pdf   # 矩形建树
├── 07_motif/
│   ├── meme_location.txt                   # MEME motif 坐标表
│   └── Tree_Domain_Motif_GeneStructure.pdf # 树 + 结构域 + motif + 外显子复合图
├── 08_jcvi_kaks/
│   └── 08_jcvi_kaks.pdf                    # Ka/Ks 分布 + 共线性 block
├── 09_circos/
│   └── circos_{species}.pdf                # 每个物种一张 Circos
├── 10_promoter/
│   └── promoter_elements.pdf               # 启动子顺式元件热图
└── 11_ppi/
    └── PPI_network.pdf                     # PPI 网络图

work/                                       # 中间文件（.aln / .nwk / .blast / fimo 等）
logs/                                       # 每条 rule 的日志
```

每个输出文件的含义见 [docs/output.md](docs/output.md)。

---

## 分步/按需运行

Snakemake 支持 `--until` 只跑到某一步、`--forcerun` 强制重跑、`--cores` 控并发：

```bash
# 只跑数据鉴定部分（Step 1-4）
snakemake --configfile config/my_config.yaml --until step04_pfam_filter -j 8 --use-conda

# 跑到建树为止
snakemake --configfile config/my_config.yaml --until step06_tree_plot -j 4 --use-conda

# 只跑 Step 10（前置步骤已完成）
snakemake --configfile config/my_config.yaml --until step10_promoter -j 2 --use-conda

# 强制重跑 Step 7（改了 MEME 参数之后）
snakemake --configfile config/my_config.yaml --forcerun step07_meme -j 4 --use-conda

# 看 DAG 图
snakemake --configfile config/my_config.yaml --dag | dot -Tpdf > dag.pdf
```

---

## ARM Mac 特别说明（Step 10 Docker）

bioconda 在 `osx-arm64`（Apple Silicon）下**没有 MEME 5.1.0**，而默认 JASPAR 束是用 5.1.0 打的。所以 M 系列 Mac 上跑 Step 10 需要切换到 `linux/amd64` 容器：

```bash
# 首次运行会自动构建镜像（~200 MB，~3 min）
scripts/run-in-linux.sh fimo --version

# 在容器里跑 Step 10
scripts/run-in-linux.sh snakemake --configfile config/my_config.yaml \
    --until step10_promoter -j 4 --use-conda

# 进交互 shell
scripts/run-in-linux.sh bash
```

镜像定义：`docker/Dockerfile.meme510`（基于 `mambaorg/micromamba:1.5-jammy`，pin `meme=5.1.0`）。

Intel Mac 和 Linux 可以直接用 conda 环境，不需要 Docker。

---

## 测试与开发

```bash
# 全量单元测试
pytest tests/ -v
# 389 passed, 1 skipped, 97% coverage

# 只跑某一个模块的测试
pytest tests/test_parse_meme_output.py -v

# dry-run 校验 Snakefile
snakemake --configfile config/default_config.yaml -n -p
```

加新分析步骤的流程：

1. `scripts/` 加 Python 工具（带类型注解和 argparse CLI）
2. `tests/test_*.py` 加单元测试
3. `R/` 加可视化脚本（optparse 接口）
4. `Snakefile` 加 rule，更新 `rule all:` 的 input
5. `config/default_config.yaml` 加对应参数块
6. `snakemake -n` 验证 DAG → 实际运行

完整贡献指南见 [docs/development.md](docs/development.md)。

---

## 文档索引

| 文档 | 内容 |
|------|------|
| [docs/tutorial.md](docs/tutorial.md) | 带案例数据的完整 step-by-step 教程 |
| [docs/configuration.md](docs/configuration.md) | `config/default_config.yaml` 每个字段的参考 |
| [docs/output.md](docs/output.md) | 每个输出文件的含义和解读方法 |
| [docs/troubleshooting.md](docs/troubleshooting.md) | 常见错误（安装/运行/外部工具） |
| [docs/development.md](docs/development.md) | 贡献者指南（代码规范、测试、新步骤） |
| [CLAUDE.md](CLAUDE.md) | 架构速览（给 Claude Code / 开发者看的） |

---

## 引用

```
Jiawang1209. GeneFamilyFlow: Automated multi-species gene family analysis pipeline.
GitHub repository, https://github.com/Jiawang1209/GeneFamilyFlow
```

## License

MIT — see [LICENSE](LICENSE).
