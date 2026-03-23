# GeneFamilyFlow

## Step0: 环境配置、编码原则和技术栈

环境配置：

- 环境需要单独的使用conda创建一个环境
- 项目会涉及到Python和R
- Python部分主要承担流程控制、R部分则主要用于生物信息学统计分析与结果展示
- 既可以实现命令行去执行，也可以分步骤执行每一步骤，最后形成一个Rmarkdown的结果报表

编码原则：

- 牢记奥卡姆剃刀原则
- 遵循第一性原理
- 避免防御性编程



----

我现在想开发一个多物种基因家族的分析流程，主要包含以下分析内容：

## Step1: 基因组数据的准备和清洗

- 我目前准备的基因组数据，目前来源于两部分，第一部分是`Phytozome`，第二部分是`1000 species`，我需要三个数据，一个是基因组序列的数据，一个是gff或者gtf文件，然后最后一个是protein.fasta。不过值得注意的是`gff`和`gtf`文件是必备的。如果能有基因组序列的话，那就更好了
- 在上面这一步数据准备i好之后，我们接下来要进行的是数据清洗，首先我们要看看基因组数据是否已经注释到了染色体，或者scaffold水平，我们可以肯定的是contig水平因该是不行的。其次我们需要基于gtf/gff3文件，以及protein.fasta文件，筛选出最长转录本，或者代表性转录本。



## Step2: hmmsearch

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/2.hmmsearch
conda activate genefamily

# 首先下载 PF00854 结构域
# https://www.ebi.ac.uk/interpro/entry/pfam/PF00854/curation/

gunzip PF00854.hmm.gz

# cp /data/Liu/LiuYue/Evolution/Evolution_LWX_NPF_2024/2.hmmsearch/PF00854.hmm ./

# 使用HMM模型搜索pep文件
ls ../1.database/*pep.fasta|while read id
do 
id1=$(basename $id)
echo "hmmsearch --cut_tc --domtblout ${id1%%.*}.PF00854.domtblout -o ${id1%%.*}.PF00854.hmmout ./PF00854.hmm ../1.database/${id1}"
done > command_hmmsearch.sh
nohup bash command_hmmsearch.sh &

# 筛选数据 基于 第七列的E-value 1e-20
ls *domtblout|while read id
do
awk '$7<1e-10 && $1 !~ /^#/ {print $0}' $id > ${id}.filter
done

# 提取序列
ls *filter|while read id
do
awk '{print $1}' $id |sort -u > ${id}.id_1st
done

cat *id_1st|sort|uniq  > 1st_id
wc -l 1st_id # 226
cat ../1.database/*pep.fasta > species_10.pep.fa
# 提取序列
seqkit grep -r -f 1st_id species_10.pep.fa -o 1st_id.fa

################
# clustalw 进行多序列比对，构建双子叶植物特定的NB-ARC结构域的HMM模型
# clustalw -> 1. Sequence Input From Disc -> 1st_id.fa # 输入数据 
# ->2. Multiple Alignments -> 1.  Do complete multiple alignment now Slow/Accurate --> 1st_id.aln -> 1st_id.dnd ->X
################
# 最终生成 1st_id.aln 文件和 1st_id.dnd 文件

nohup clustalw -infile=1st_id.fa -output=clustal -type=PROTEIN -outfile=1st_id.aln &

# hmmbuild构建新的模型
hmmbuild new_NPF.hmm 1st_id.aln

# 使用新的HMM模型重新检索pep文件
ls ../1.database/*pep.fasta|while read id
do 
id1=$(basename $id)
echo "hmmsearch --domtblout ${id1%%.*}.new_NPF.domtblout -o ${id1%%.*}.new_NPF.hmmout ./new_NPF.hmm ../1.database/${id1}"
done > command_hmmsearch_2st.sh
nohup bash command_hmmsearch_2st.sh &

# 获取hmmsearch的数据

ls *new_NPF.domtblout|while read id
do
awk '$7<1e-10 && $1 !~ /^#/ {print $0}' $id > ${id}.filter_2st
done

ls *filter_2st|while read id
do
awk '{print $1}' $id |sort -u > ${id}.id_2st
done

cat *id_2st > 2st_id
wc -l 2st_id # 239
```



## Step3:blast

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/3.blast
# 基于 AT 寻找 NBS-LRR的序列，进行blastp
cp /data/Liu/LiuYue/Evolution/Gramineae_Evolution/MTP_21/3_blast/all.domains.txt ./
cp /data/Liu/LiuYue/Evolution/Gramineae_Evolution/MTP_21/3_blast/AT.clean.pep.fasta ./

grep 'PF00854' all.domains.txt|awk -F '.' '{print $1}'|sort|uniq > PF00854.TAIR.ID

# 53

# 获取序列
seqkit grep -r -f PF00854.TAIR.ID AT.clean.pep.fasta -o PF00854.TAIR.ID.fa

# blastp
makeblastdb -in PF00854.TAIR.ID.fa -dbtype prot -title -parse_seqids

ls ../1.database/*pep.fasta|while read id
do 
id1=$(basename $id)
nohup blastp -query ../1.database/${id1%%.*}.pep.fasta -db PF00854.TAIR.ID.fa -outfmt '6 std qlen slen' -out ./${id1%%.*}.blast -evalue 1e-10 -num_threads 10 -num_alignments 10 & 
done

# 筛选
cat *blast > species_10.blast
# awk '$3 > 30 && ($8-$7+1) / $13 > 0.5 {print $1}' species_10.blast|sort|uniq > species_10.blast.id
awk '{print $1}' species_10.blast|sort|uniq > species_10.blast.id


wc -l species_10.blast.id # 240
```



## Step4:identification

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/4.identification

# 把两个结果进行交集
cat ../2.hmmsearch/2st_id ../3.blast/species_10.blast.id | sort |uniq -c |awk '$1 == 2{print $2}' > inter.ID

wc -l inter.ID # 238

cat ../2.hmmsearch/2st_id ../3.blast/species_10.blast.id|uniq -c > union.ID

wc -l union.ID # 479

seqkit grep -r -f inter.ID ../2.hmmsearch/species_10.pep.fa -o inter.ID.fa

# pfam 鉴定
nohup hmmscan -o pfam.out.txt --tblout pfam.out.tbl --cpu 30 --noali -E 1e-5 ~/software/hmmer-3.1b1/Pfam-A.hmm ./inter.ID.fa &

grep 'PF00854' pfam.out.tbl|awk '{print $3}'|sort|uniq|wc -l # 226

grep 'PF00854' pfam.out.tbl|awk '{print $3}'|sort|uniq > pfam.ID

seqkit grep -r -f pfam.ID ../2.hmmsearch/species_10.pep.fa -o pfam.ID.fa


# pfam_scam.pl 鉴定
conda activate rna
nohup pfam_scan.pl -fasta inter.ID.fa -dir /data/Liu/LiuYue/Evolution/Pfam/ -cpu 30 -out Pfam_scan.out > Pfam_scan.log 2>&1 &


grep 'PF00854' Pfam_scan.out|awk '{print $1}'|sort|uniq|wc -l # 1109
grep 'PF00854' Pfam_scan.out|awk '{print $1}'|sort|uniq > pfam_scan.id


seqkit grep -r -f pfam_scan.id ../2.hmmsearch/species_10.pep.fa -o identify.ID.fa

# 总个数：1109

# clean 
awk '{print $1}' identify.ID.fa > identify.ID.clean.fa
```





## Step5:genefamily_info

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/5.genefamily_info

cp ../4.identification/identify.ID.fa ./

awk '{print $1}' identify.ID.fa > original.ID.fa

# clean 
awk '{print $1}' identify.ID.fa > identify.ID.clean.fa

# 准备 clean ID
awk '{print $1}' identify.ID.fa |grep '>'|sed 's/>//g'|sort|uniq > original.ID

cp ../../Evolution_Jiangxi_Liu_2024/5.genefamily_info/clean_ID.R ./

# conda activate DataScience
Rscript clean_ID.R original.ID
cat ../1.database/*gff* > species_9.gff3


awk '$3=="gene" {print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}' species_9.gff3|awk -F ';' '{print $1}' | sed 's/ID=//g'|awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > species_10.bed
```



## Step6:tree

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/6.tree

# 多序列比对
nohup muscle -in ../4.identification/identify.ID.clean.fa   -out  identify.ID.muscle 2> identify.ID.muscle.log &

# iqtree
nohup iqtree -s identify.ID.muscle -m MFP -bb 1000 -bnni  -nt AUTO  -cmax 15  -redo  &

grep '>So' ../4.identification/identify.ID.clean.fa|sed 's/>//g' > Sb.NPF.ID
seqkit grep -r -f Sb.NPF.ID ../4.identification/identify.ID.clean.fa -o Sb.NPF.ID.fa


# 多序列比对
nohup muscle -in Sb.NPF.ID.fa   -out  Sb.NPF.ID.muscle 2> Sb.NPF.ID.muscle.log &

# iqtree
nohup iqtree -s Sb.NPF.ID.muscle -m MFP -bb 1000 -bnni  -nt AUTO  -cmax 15  -redo  &
```



## Step7:motif_genestructure

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/7.motif_genestructure
cp ../6.tree/Sb.NPF.ID.fa ./
nohup meme Sb.NPF.ID.fa  -mod anr  -protein -nmotifs 10 -minw 10 -maxw 50 &
```





## Step8:JCVI

```shell
# Osativa
# Athaliana
# Sbicolor

cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/8.collinearity

# 准备文件
cp ../1.database/*gff3 ./
cp ../1.database/*pep.fasta ./

# 转换 GFF 文件到 BED文件 并且重新命名
ls *gff3|while read id
do
nohup python -m jcvi.formats.gff bed --type=mRNA --key=ID ${id} -o ${id%%.*}.bed &
done

# bed 文件 去除重复
ls *bed|while read id
do
nohup python -m jcvi.formats.bed uniq $id &
done

# 然后我们对 uniq.bed 文件再次进行清洗
ls *uniq.bed|while read id
do 
echo ${id}
head -n 10 $id
done

# ls *uniq.bed|while read id; do  echo ${id}; head -n 10 $id; done

# Athaliana.uniq.bed 不需要修改
# Osativa.uniq.bed 需要修改， 修改 Osativa.uniq.bed 即可
perl -p -i -e 's/.MSUv7.0//g' Osativa.uniq.bed

# Sbicolor.uniq.bed 修改 Sbicolor.uniq.bed 即可
perl -p -i -e 's/.v3.1//g' Sbicolor.uniq.bed
# 也要同时修改一下 Sbicolor.pep.fasta
perl -p -i -e 's/\.p.*//g' Sbicolor.pep.fasta

# 准备 pep 序列
ls *uniq.bed|while read id
do
id1=${id%%.*}
seqkit grep -f <(cut -f 4 ${id1}.uniq.bed ) ${id1}.pep.fasta | seqkit seq -i > ${id1}.pep
done

# 开始分析数据
mkdir -p out && cd out
ln -s ../*pep ./

ls ../*uniq.bed|while read id
do
id1=$(basename ${id})
id2=${id1%%.*}
ln -s ${id} ${id2}.bed
done

# 然后开始共线性分析
# 共线性分析
python -m jcvi.compara.catalog ortholog --dbtype prot --notex  --no_strip_names Sbicolor Athaliana 
python -m jcvi.compara.catalog ortholog --dbtype prot --notex  --no_strip_names Athaliana Osativa
python -m jcvi.compara.catalog ortholog --dbtype prot --notex  --no_strip_names Osativa Sbicolor
```

> 创建配置文件 `.simple`文件

```shell
ls *anchors| grep -v 'lifted'| while read id
do
python -m jcvi.compara.synteny screen --minspan=30 --simple ${id} ${id}.new
done
```

> 创建配置文件：`seqids`文件

```shell
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10
Chr1,Chr2,Chr3,Chr4,Chr5
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12
```

> 创建配置文件：`layout`文件

```shell
# y, xstart, xend, rotation, color, label, va,  bed
.8,      0.2,    .8,      0,     #dd3497, Sbicolor, top, Sbicolor.bed
.5,      0.025,    0.625,      120,      #4eb3d3, Athaliana, top, Athaliana.bed
.5,     0.375,    0.975,     60,      #807dba, Osativa, top, Osativa.bed
# edges
e, 0, 1, Sbicolor.Athaliana.anchors.simple
e, 1, 2, Athaliana.Osativa.anchors.simple
e, 0, 2, Osativa.Sbicolor.anchors.simple
```

> 添加颜色

```shell
# apple id
grep '>' ../../4.identification/identify.ID.fa| awk '{print $1}'| sed 's/>//g'| sed 's/.p$//g'| sed -E 's/.cds.*//g' | sed -E 's/_P0.*//g' | sed -E 's/.P0.*$//g' > select.ID

wc -l select.ID # 228


# 构造数据
ls *anchors.new| while read id
do
grep -Ff select.ID -w $id > ${id%.*}.color
done

# 构造数据
# ls *color|while read id
# do
# awk '{print "g*"$1"\t"$1"\t"$2"\t"$2"\t"$3"\t""+"}' $id > ${id}2
# done
awk '{print "r*"$1"\t"$1"\t"$2"\t"$2"\t"$3"\t""+"}' Athaliana.Osativa.anchors.color > Athaliana.Osativa.anchors.color2
awk '{print "g*"$1"\t"$1"\t"$2"\t"$2"\t"$3"\t""+"}' Osativa.Sbicolor.anchors.color > Osativa.Sbicolor.anchors.color2
awk '{print "b*"$1"\t"$1"\t"$2"\t"$2"\t"$3"\t""+"}' Sbicolor.Athaliana.anchors.color > Sbicolor.Athaliana.anchors.color2



# 合并数据
ls *anchors.simple|while read id
do
cat ${id} ${id%.*}.color2 > ${id%.*}.simple2
mv ${id%.*}.simple2 ${id%.*}.simple
done
```

> 最后绘图

```shell
python -m jcvi.graphics.karyotype seqids layout --notex --figsize=14x12 --chrstyle=roundrect
```



### Step8.1: KaKs

> KaKs 计算和统计

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/8.collinearity/out

# 共线性的数据
cat *color2 > jcvi_kaks_all.file

wc -l jcvi_kaks_all.file 
# 48 jcvi_kaks_all.file

# 准备 gene-pair
awk -F '*' '{print $2}' jcvi_kaks_all.file|awk '{print $1"\t"$3}' > KaKs_Gene_Pair

wc -l KaKs_Gene_Pair
# 48 KaKs_Gene_Pair


# 准备ID
awk -F '*' '{print $2}' jcvi_kaks_all.file|awk '{print $1}'|sort|uniq > Kaks_Gene_ID
awk -F '*' '{print $2}' jcvi_kaks_all.file|awk '{print $3}'|sort|uniq >> Kaks_Gene_ID
cat Kaks_Gene_ID|sort|uniq > Kaks_Gene_ID2
mv Kaks_Gene_ID2 Kaks_Gene_ID
wc -l Kaks_Gene_ID
# 78 Kaks_Gene_ID

# cds fasta 需要额外清洗数据
# Athaliana
awk '{print $1}' ../../1.database/Athaliana.cds.fasta > Athaliana.cds.fasta

# Osativa
awk '{print $1}' ../../1.database/Osativa.cds.fasta > Osativa.cds.fasta

# Sbicolor
awk '{print $1}' ../../1.database/Sbicolor.cds.fasta | perl -p -e 's/\.p.*//g' > Sbicolor.cds.fasta

# 合并数据
cat *cds.fasta|awk '{print $1}' > species_10.cds.fa
seqkit grep -f Kaks_Gene_ID  species_10.cds.fa -o Kaks_Gene_ID.cds.fasta

# check
grep '>' Kaks_Gene_ID.cds.fasta | wc -l # 78

```



## Step9:mcscanx

```shell
#### 最终保留如下物种 ####
# Osativa
# Sbicolor
# Athaliana

/home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/9.mcscanx

# 准备数据 pep.fa
ln -s ../8.collinearity/Sbicolor.pep ./
ln -s ../8.collinearity/Athaliana.pep ./
ln -s ../8.collinearity/Osativa.pep ./

# blastp 
nohup makeblastdb -in Sbicolor.pep -dbtype prot -title -parse_seqids &
nohup makeblastdb -in Athaliana.pep -dbtype prot -title -parse_seqids &
nohup makeblastdb -in Osativa.pep -dbtype prot -title -parse_seqids &

nohup blastp -query Sbicolor.pep -db Sbicolor.pep -outfmt 6 -out ./Sbicolor.blast -num_threads 15 -num_alignments 5 -evalue 1e-10  & 
nohup blastp -query Athaliana.pep -db Athaliana.pep -outfmt 6 -out ./Athaliana.blast -num_threads 15 -num_alignments 5 -evalue 1e-10  & 
nohup blastp -query Osativa.pep -db Osativa.pep -outfmt 6 -out ./Osativa.blast -num_threads 15 -num_alignments 5 -evalue 1e-10  & 

cp ../8.collinearity/Athaliana.bed ./
cp ../8.collinearity/Sbicolor.bed ./
cp ../8.collinearity/Osativa.bed ./

awk '{print $1"\t"$4"\t"$2"\t"$3}' ../8.collinearity/Sbicolor.uniq.bed > Sbicolor.gff
awk '{print $1"\t"$4"\t"$2"\t"$3}' ../8.collinearity/Athaliana.uniq.bed > Athaliana.gff
awk '{print $1"\t"$4"\t"$2"\t"$3}' ../8.collinearity/Osativa.uniq.bed > Osativa.gff

# MCScanX
~/software/MCScanX/MCScanX ./Sbicolor
~/software/MCScanX/MCScanX ./Athaliana
~/software/MCScanX/MCScanX ./Osativa

~/software/MCScanX/duplicate_gene_classifier ./Sbicolor
~/software/MCScanX/duplicate_gene_classifier ./Athaliana
~/software/MCScanX/duplicate_gene_classifier ./Osativa

# 准备染色体长度文件
cp /data/Liu/LiuYue/Evolution/Gramineae_Evolution/1.Gramineae_database/genome_file/genome_length/Athaliana*fai ./
cp /data/Liu/LiuYue/Evolution/Gramineae_Evolution/1.Gramineae_database/genome_file/genome_length/Osativa*fai ./
cp /data/Liu/LiuYue/Evolution/Gramineae_Evolution/1.Gramineae_database/genome_file/genome_length/Sbicolor*fai ./

# 准备ID
grep '>' ../5.genefamily_info/identify.ID.clean.fa|grep 'So'|sed 's/>//g' > Sbicolor.NPF.id
grep '>' ../5.genefamily_info/identify.ID.clean.fa|grep 'AT'|sed 's/>//g' > AT.NPF.id
grep '>' ../5.genefamily_info/identify.ID.clean.fa|grep 'LOC'|sed 's/>//g' > Os.NPF.id

# 数据已经都准备好了，然后开始画图
```

### Step9.1:KaKs

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/9.mcscanx

# 准备cds 数据
cp ../8.collinearity/out/species_10.cds.fa ./

# 准备ID
awk '{print $1}'  mcscanx_Ka_Ks.pair2.txt > kaks_1.ID
awk '{print $2}'  mcscanx_Ka_Ks.pair2.txt > kaks_2.ID
cat kaks_1.ID kaks_2.ID |sort|uniq > kaks.ID

rm kaks_1.ID kaks_2.ID

seqkit grep -f kaks.ID species_10.cds.fa -o kaks.ID.cds.fasta

grep '>' kaks.ID.cds.fasta |wc -l # 184

```





## Step10:promoter

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/10.promoter
# 准备ID文件
cp ../6.tree/Sb.NPF.ID ./
perl -p -i -e 's/\.\d\.p//g' Sb.NPF.ID


# 准备GFF3文件
cat ../1.database/*gff3 > species_10.gff3

# 清洗数据
awk '$3=="gene" {print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}' species_10.gff3|awk -F ';' '{print $1}' | sed 's/ID=//g'|awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > species_10.bed

# 提取候选家族成员的bed
grep -Ff Sb.NPF.ID -w species_10.bed > original.ID.clean.bed


# 提取序列
# Sb
seqkit subseq --bed original.ID.clean.bed --up-stream 2000 --only-flank /data/Liu/LiuYue/Evolution/Gramineae_Evolution/1.Gramineae_database/genome_file/Sbicolor_454_v3.0.1.fa.gz > Sb.promoter.fasta
```



## Step11:ppi

```shell
cd /home/yhc-z4/workspace/Evolution/Evolution_NPF_Sb/11.ppi

# 准备好数据
cp ../5.genefamily_info/identify.ID.clean.fa ./

# 分别提取出ID 用于进行序列比对
grep '>AT'  identify.ID.clean.fa|sed 's/>//g' > AT.pep.ID
grep '>LOC' identify.ID.clean.fa|sed 's/>//g' > Os.pep.ID
grep '>So' identify.ID.clean.fa|sed 's/>//g' > Sbicolor.pep.ID

# 获取序列
ls *pep.ID|while read id
do
seqkit grep -f <(cut -f 1 ${id}) identify.ID.clean.fa -o ${id%%.*}.GF.pep.fasta
done

# check
ls *pep.ID|while read id
do 
echo ${id}
wc -l ${id}
grep '>' ${id%%.*}.GF.pep.fasta|wc -l
done

# 准备蛋白质序列
ln -s ../1.database/*pep.fasta ./


ls *pep.fasta|grep -v 'GF'|while read id
do
makeblastdb -in ${id} -dbtype prot -title -parse_seqids
done

# 所有的家族基因蛋白 都要和 Athaliana.pep.fasta 去做比对
ls *pep.fasta|grep 'GF'|while read id
do
nohup blastp -query ${id} -db Athaliana.pep.fasta -outfmt '6 std qlen slen' -out ./${id%%.*}_Athaliana.blast -evalue 1e-10 -num_threads 10 -num_alignments 3 &
done

# Athaliana.pep.fasta 又要和所有的其他物种的蛋白去做比对
ls *pep.fasta|grep -v 'GF'|while read id
do
nohup blastp -query Athaliana.pep.fasta -db ${id} -outfmt '6 std qlen slen' -out ./Athaliana_${id%%.*}.second.blast -evalue 1e-10 -num_threads 10 -num_alignments 3 &
done
```



>从`PPI`的结果里面倒出`node_annotation_tmp.txt`

```shell
awk '{print $1}'  ../2.hmmsearch/species_10.pep.fa  > species_10.pep.fa

grep '>' species_10.pep.fa|sed 's/>//g' > species_10.pep.ID

grep -Ff PPI.geneid.txt -w species_10.pep.ID > all.ID.pep.select

seqkit grep -f all.ID.pep.select species_10.pep.fa -o all.ID.pep.select.pep.fasta

conda activate rna
nohup pfam_scan.pl -fasta all.ID.pep.select.pep.fasta -dir /data/Liu/LiuYue/Evolution/Pfam/ -cpu 30 -out all.ID.pep.select.Pfam_scan.out > all.ID.pep.select.Pfam_scan.log 2>&1 &
```





























