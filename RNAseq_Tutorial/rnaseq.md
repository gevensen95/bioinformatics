Introduction to RNAseq Analysis
================
K. Garrett Evensen
5/10/2022

## Data Preparatation

Analyzing RNAseq data may seem daunting, but after this tutorial you
will have a good understanding of the basics of how RNAseq data is
handled and analyzed. The data for this tutorial comes from two
peer-reviewed papers (<https://doi.org/10.1016/j.cbd.2021.100952>,
<https://doi.org/10.1016/j.margen.2021.100865>) on the blue mussel,
Mytilus edulis. One is focused on identification of putatitve sex
diffferentiating genes using phylotranscriptomics and the other is on
evaluating expression differences in deep and shallow water growth.
Here, we will re-create the differential gene expression analysis
section of the former paper.

The first input file for this tutorial is a counts matrix. This is file
with samples as columns and genes as headers. If you open the file, you
will see a bunch of number for each gene and sample. These are called
counts, which we will use to identify differentially expressed genes
(DEG). I made this file by aligning the trimmed and quality-filtered
read files for each sample to a transcriptome I assembled. I won’t be
covering these steps in this tutorial, but I will paste some example
code for trimming reads, assembling a transcriptome, and calculating
gene counts that you can use as a starting point with your own RNAseq
data. COUDBLE CHECK PAIRED OPTION FOR TRINITY

trim_galore –stringency 3 -e 0 –paired forward.fastq.gz reverse.fastq.gz
–fastqc

Trinity –seqType fq –max_memory 50G –SS_lib_type F –normalize_reads –CPU
20 –paired ./merged_trimmed_forward.fq.gz –min_contig_length 400 –output
./Trinity_out

rsem-prepare-reference –bowtie2 -p 4 Trinity.fasta Trinity
rsem-calculate-expression –bowtie2 -p 7 –paired-end
Sample1_forward.fastq.gz Sample1_reverse.fastq.gz Trinity Hg7d_3

In this tutorial, we will use limma to determine DEG. There are many
different methods for DEG analysis, such as DESeq2 and EBSeq.

If you don’t have limma or edgeR installed, run this code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.18), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'limma'

    ## Old packages: 'aplot', 'ellipse', 'geojsonsf', 'openssl', 'usethis'

``` r
BiocManager::install("edgeR")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.18), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'edgeR'

    ## Old packages: 'aplot', 'ellipse', 'geojsonsf', 'openssl', 'usethis'

``` r
require(limma)
```

    ## Loading required package: limma

``` r
require(edgeR)
```

    ## Loading required package: edgeR

``` r
require(tidyverse)
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
counts <- data.matrix(read.csv("counts.csv", row.names = 1, header=T))
annot <- read.csv("annotation.csv", row.names = 1, header = T)
```

Now that we have imported our data, we need to make a design matrix.
This tells limma what our experiment is. THIS NEEDS TO SOUND BETTER. We
have 12 females and 12 males in our dataset. Luckily, they are in order,
with females first. In our design matrix, we have set females as the
intercept. If you have a control group, then they would be set as the
intercept. If you use a contrast matrix (see below), then you would use
“design \<- model.matrix(\~0+conditions)” instead. This ensures that no
group is set as the intercept.

``` r
f <- rep("F", 12)
m <- rep("M", 12)
conditions = factor(c(f,m))

design <- model.matrix(~conditions) 
colnames(design) <- levels(conditions)
design
```

    ##    F M
    ## 1  1 0
    ## 2  1 0
    ## 3  1 0
    ## 4  1 0
    ## 5  1 0
    ## 6  1 0
    ## 7  1 0
    ## 8  1 0
    ## 9  1 0
    ## 10 1 0
    ## 11 1 0
    ## 12 1 0
    ## 13 1 1
    ## 14 1 1
    ## 15 1 1
    ## 16 1 1
    ## 17 1 1
    ## 18 1 1
    ## 19 1 1
    ## 20 1 1
    ## 21 1 1
    ## 22 1 1
    ## 23 1 1
    ## 24 1 1
    ## attr(,"assign")
    ## [1] 0 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$conditions
    ## [1] "contr.treatment"

Next, we will filter out any genes with low expression across all our
samples. In this case, NEED TEXT FROM PAPER. Then we will normalize our
counts using the TMM method and calculate the log counts per million
(logCPM) of our gene counts.

``` r
dge <- DGEList(counts = counts, genes = annot) 
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM") 
logCPM <- cpm(dge, log=TRUE, prior.count = 3)
```

So, how many genes were actually kept?

``` r
length(counts[,1])
```

    ## [1] 98462

``` r
length(logCPM[,1])
```

    ## [1] 23958

If we compare those two numbers, we started with 98,462 genes, and
74,504 genes were removed due to low expression values, leaving 23,958
genes. For reference, there are 65,625 genes in the Mytilus
galloprovincialis genome. The number of genes in the Mytilus edulis
genome is not presently known as the genome is still being annotated.
Here we are using an assembled transcriptome rather than a genome
because Mytilus sp. genome’s are extremely hemizgous and the M. edulis
genome is not annotated, yet.

Sometimes, the last step before we determine DEG is creating a contrast
matrix. This is used to determine what comparisons to make and if genes
are up in a particular condition. The code below is for if we had a
third condition. Because we only have two conditions, we won’t use that
for any analysis.

contrast.matrix \<- makeContrasts(F - M, H - F, M - H, levels = design)

# Analysis

Now we can jump into the analysis steps! Limma becomes much more
powerful when voom is added.

``` r
v <- voom(dge, design, plot=TRUE)
```

![](rnaseq_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
fit <- lmFit(v, design)
fit.eBayes = eBayes(fit, trend = TRUE)

summary <- summary(decideTests(fit.eBayes))
summary
```

    ##            F     M
    ## Down    6267  4732
    ## NotSig  6060 13713
    ## Up     11631  5513

If we used a contrast matrix, we would add the code bellow to determine
which genes are differentially expressed and how many for each
comparison.

fit2 \<- contrasts.fit(fit, contrast.matrix) fit2.eBayes = eBayes(fit2,
trend = TRUE) summary2 \<- summary(decideTests(fit2.eBayes))

Because we set females (F) as the intercept, we will look only at the
number of DEG in the male (M) category. Overall, 4,732 genes are
down-regulated and 5,513 genes are up-regulated in males.

We’ve done the analysis, but now we need to create a data frame that has
all the genes, their fold changes, and their adjusted p values.

``` r
male <- topTreat(fit.eBayes, coef = "M", n = Inf, sort = "p", adjust.method = "BH")
```

Now you have a list of all the genes, their annotations, their log2 fold
change, BH adjusted p vales, and a few other important pieces of
information.

If you only want the DEG, then you can add p.value to topTreat.

``` r
male2 <- topTreat(fit.eBayes, coef = "M", n = Inf, sort = "p", adjust.method = "BH", p.value = 0.05)
```

Next, we will visualize the number of DEG in males.

``` r
threshold_m <- male$adj.P.Val < 0.05 
ggplot(male) +
  geom_point(aes(x=logFC, y=-log10(adj.P.Val), color=threshold_m)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,30)) +
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.25), face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(face="bold",vjust = .5),
        axis.text.y = element_text(face = "bold")) +
  scale_color_discrete(name = "DEG", labels = c("Yes", "No"))
```

![](rnaseq_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Heatmap

The next step is to create a heatmap to see if there are any overlapping
samples.

``` r
require(gplots)
```

    ## Loading required package: gplots

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
sig_level <- 0.05
male_sig <- rownames(male %>%
  filter(adj.P.Val < sig_level))
logCPM <- as.data.frame(logCPM)

heatmap_data <- logCPM %>%
  filter(rownames(logCPM) %in% male_sig)
heatmap_data <- data.matrix(heatmap_data)
heatmap.2(heatmap_data, scale = "row", trace = "none", cexCol = 1, lhei = c(1,3),
          labRow = FALSE)
```

![](rnaseq_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# PCA

Finally, we will use PCA on the DEG to see if the male and female
samples can be separated. We will use the logCPM object we made earlier
to do this.

``` r
logCPM_r <- as.data.frame(t(logCPM)) #rotates it so sites are rows and genes are columns
sex <- c(rep("F",12), rep("M", 12))
logCPM_r$sex <- sex
pca <- prcomp(logCPM_r[,-23959]) #PCA on only gene counts
summary(pca)
```

    ## Importance of components:
    ##                             PC1     PC2      PC3      PC4      PC5      PC6
    ## Standard deviation     138.6205 68.4425 48.50624 44.24086 32.90876 32.37913
    ## Proportion of Variance   0.4554  0.1110  0.05576  0.04638  0.02566  0.02484
    ## Cumulative Proportion    0.4554  0.5664  0.62213  0.66851  0.69418  0.71902
    ##                             PC7      PC8      PC9     PC10    PC11     PC12
    ## Standard deviation     30.16649 29.18640 28.58260 28.15361 27.5579 27.21147
    ## Proportion of Variance  0.02157  0.02019  0.01936  0.01878  0.0180  0.01755
    ## Cumulative Proportion   0.74059  0.76078  0.78014  0.79892  0.8169  0.83446
    ##                            PC13     PC14     PC15     PC16     PC17     PC18
    ## Standard deviation     27.14075 26.68563 26.27455 25.74863 25.31543 25.14255
    ## Proportion of Variance  0.01746  0.01688  0.01636  0.01571  0.01519  0.01498
    ## Cumulative Proportion   0.85192  0.86880  0.88516  0.90087  0.91605  0.93103
    ##                            PC19     PC20     PC21     PC22     PC23     PC24
    ## Standard deviation     24.86872 24.25935 23.98555 23.82321 23.67325 1.67e-13
    ## Proportion of Variance  0.01466  0.01395  0.01363  0.01345  0.01328 0.00e+00
    ## Cumulative Proportion   0.94569  0.95964  0.97327  0.98672  1.00000 1.00e+00

``` r
plot(pca)
```

![](rnaseq_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

As you can see, most of the variance is in PC1 (45.54%)

``` r
require(factoextra)
```

    ## Loading required package: factoextra

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
fviz_pca_ind(pca, axes = c(1,2), 
             geom.ind = "point",
             col.ind = logCPM_r$sex,
             legend.title = "Sex",
             xlab = "PC1 (45.54%)", ylab = "PC2 (11.1%)",
             pointsize = 5,
             title = "", theme_classic()
)
```

![](rnaseq_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Based on this PCA plot using PC1 and PC2, there is a clear separation
between male and female samples.

# GO & KEGG Enrichment

An important aspect of RNAseq analysis is understanding which pathways
are affected by a particular treatment. We will investigate this by
using both GO and KEGG enrichment. I have already collected GO terms and
KO ID’s for our genes using the Trinotate pipeline and BLAST Koala. We
will start with GO enrichment using the topGO package.

``` r
BiocManager::install("topGO")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.18), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'topGO'

    ## Old packages: 'aplot', 'ellipse', 'geojsonsf', 'openssl', 'usethis'

``` r
require(topGO)
```

    ## Loading required package: topGO

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: graph

    ## 
    ## Attaching package: 'graph'

    ## The following object is masked from 'package:stringr':
    ## 
    ##     boundary

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: GO.db

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:gplots':
    ## 
    ##     space

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    ## 
    ## groupGOTerms:    GOBPTerm, GOMFTerm, GOCCTerm environments built.

    ## 
    ## Attaching package: 'topGO'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     members

Now that we have topGO installed and loaded, we will reformat our data
to investigate enrichmed GO terms among the upregulated genes.

``` r
up_genes <- male %>%
  filter(logFC > 0)
up_genes_names <- rownames(up_genes)
up_gene_p <- up_genes$adj.P.Val
names(up_gene_p) <- up_genes_names
```

Next, we can do GO enrichment!

``` r
geneID2GO <- readMappings("gene_go.csv", sep = ",", IDsep = ";")
DEG<-function (allScore) {
  return(allScore < 0.05)
}

male_topGO <- new("topGOdata", ontology = "BP", geneSelectionFun=DEG, allGenes = up_gene_p, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
```

    ## 
    ## Building most specific GOs .....

    ##  ( 11577 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 12392 GO terms and 27800 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 6533 genes annotated to the GO terms. )

``` r
male_KS <- runTest(male_topGO, algorithm = "classic", statistic = "ks")
```

    ## 
    ##           -- Classic Algorithm -- 
    ## 
    ##       the algorithm is scoring 12392 nontrivial nodes
    ##       parameters: 
    ##           test statistic: ks
    ##           score order: increasing

``` r
male_results <- GenTable(male_topGO, KS = male_KS,
                         orderBy = "KS", ranksOf = "classic", topNodes = 10)
require(ggpubr)
```

    ## Loading required package: ggpubr

``` r
go_table <- ggtexttable(male_results, rows = NULL)
go_table
```

![](rnaseq_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Now that we have finished GO enrichment in the upregulated genes, we
will do the same for KEGG enrichment.

``` r
BiocManager::install('clusterProfiler')
```

    ## Bioconductor version 3.15 (BiocManager 1.30.18), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'clusterProfiler'

    ## Old packages: 'aplot', 'ellipse', 'geojsonsf', 'openssl', 'usethis'

``` r
require(clusterProfiler)
```

    ## Loading required package: clusterProfiler

    ## clusterProfiler v4.4.1  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
keggs_all_genes <- read_csv('keggs_all_genes.csv') 
```

    ## Rows: 32325 Columns: 2

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): geneID, KO
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
up_genes2 <- male[,c(3,7)]
up_genes2$geneID <- rownames(up_genes2)

keggs_genes <- merge(keggs_all_genes, up_genes2, by = 'geneID')
head(keggs_genes)
```

    ##                   geneID     KO     logFC  adj.P.Val
    ## 1 TRINITY_DN100009_c0_g1 K19909 2.0132219 0.00649925
    ## 2 TRINITY_DN100010_c0_g1 K12303 0.9802509 0.03848415
    ## 3 TRINITY_DN100010_c0_g1 K12303 0.9802509 0.03848415
    ## 4 TRINITY_DN100010_c0_g1 K12303 0.9802509 0.03848415
    ## 5 TRINITY_DN100010_c0_g1 K12303 0.9802509 0.03848415
    ## 6 TRINITY_DN100010_c0_g1 K12303 0.9802509 0.03848415

Now that we have combined these two datasets, we can do KEGG enrichment
on the upregulated genes.

``` r
keggs <- keggs_genes %>%
  filter(logFC > 0)
keggs <- keggs[,c(2,4)] 

kegg_en <- keggs %>%
  filter(adj.P.Val < 0.05)
kegg_en <- kegg_en$KO
kegg_en <- as.character(kegg_en)
kk <- enrichKEGG(gene         = kegg_en,
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
```

    ## Reading KEGG annotation online:
    ## 
    ## Reading KEGG annotation online:

``` r
head(kk)
```

    ##              ID                                       Description GeneRatio
    ## ko05014 ko05014                     Amyotrophic lateral sclerosis   120/992
    ## ko05016 ko05016                                Huntington disease   100/992
    ## ko05022 ko05022 Pathways of neurodegeneration - multiple diseases   122/992
    ## ko05012 ko05012                                 Parkinson disease    87/992
    ## ko05020 ko05020                                     Prion disease    85/992
    ## ko05010 ko05010                                 Alzheimer disease    92/992
    ##           BgRatio       pvalue     p.adjust       qvalue
    ## ko05014 280/13290 5.148742e-62 1.848398e-59 1.316994e-59
    ## ko05016 229/13290 1.508470e-52 2.707704e-50 1.929254e-50
    ## ko05022 368/13290 1.733454e-48 2.074367e-46 1.477998e-46
    ## ko05012 209/13290 1.085533e-43 9.742662e-42 6.941700e-42
    ## ko05020 209/13290 9.653517e-42 6.931225e-40 4.938536e-40
    ## ko05010 300/13290 1.641568e-33 9.822049e-32 6.998264e-32
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        geneID
    ## ko05014               K03968/K03883/K10410/K05759/K11518/K00417/K00416/K03066/K03065/K02127/K10408/K04432/K00235/K02729/K08738/K03035/K04570/K03037/K05692/K03030/K03943/K00237/K07874/K20858/K00432/K11519/K05456/K12896/K03039/K03954/K03033/K11143/K00415/K17908/K03953/K12192/K07374/K04523/K03934/K02725/K17095/K10409/K03965/K23625/K03032/K03062/K03064/K23609/K03061/K04348/K23600/K23611/K06268/K12882/K09490/K14284/K08331/K12741/K12881/K03028/K05688/K13525/K08269/K10396/K07375/K03880/K03881/K03879/K05410/K17589/K22913/K10706/K04648/K03029/K14296/K14297/K09291/K17906/K04575/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K02138/K03939/K02739/K03937/K02726/K08334/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K10411/K03942/K02734/K10412/K00413/K11351/K02264/K03941/K03063/K03940/K02730/K02728/K02732
    ## ko05016                                                                                                                                                           K03968/K03883/K10410/K00417/K00416/K03066/K03065/K02127/K10408/K00235/K02729/K08738/K03035/K03037/K03030/K03943/K00237/K00432/K03039/K03954/K03033/K11143/K03129/K00415/K17908/K05870/K03953/K07374/K11831/K04958/K03934/K02725/K10409/K03965/K05863/K03006/K03032/K03062/K03064/K03061/K03014/K08331/K04646/K04560/K03028/K08269/K10396/K07375/K03880/K03881/K03879/K17589/K04498/K11824/K11644/K04648/K03029/K17906/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K02138/K03939/K02739/K03937/K02726/K08334/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K10411/K03942/K02734/K10412/K00413/K11351/K02264/K04638/K03941/K03063/K03940/K02730/K02728/K02732
    ## ko05022 K03968/K03883/K10410/K11518/K00417/K00416/K14021/K03066/K03065/K02127/K10408/K02183/K04432/K00235/K02729/K08738/K03035/K04570/K03037/K08957/K03030/K03943/K00237/K07874/K05687/K02580/K20858/K00432/K11519/K03039/K03954/K03033/K11143/K00415/K17908/K03953/K12192/K07374/K04958/K03934/K02725/K10409/K03965/K05863/K23625/K03032/K03062/K03064/K23609/K03061/K04348/K23600/K23611/K06268/K09490/K08331/K04560/K03028/K03097/K05688/K13525/K08269/K10396/K07375/K03880/K03881/K03879/K03083/K10575/K04552/K05410/K17589/K22913/K03178/K08844/K04648/K03029/K04962/K17906/K04575/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K02138/K03939/K02739/K03937/K02726/K08334/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K10411/K03942/K02734/K10412/K00413/K11351/K02264/K04638/K03941/K03063/K03940/K02730/K02728/K02732
    ## ko05012                                                                                                                                                                                                                                                      K03968/K03883/K00417/K00416/K03066/K03065/K02127/K02183/K00235/K02729/K08738/K03035/K04570/K03037/K03030/K03943/K00237/K05687/K20858/K03039/K03954/K03033/K00415/K03953/K07374/K04958/K03934/K02725/K03965/K05863/K03032/K03062/K03064/K03061/K09490/K03028/K05688/K10396/K07375/K03880/K03881/K03879/K10575/K04552/K04345/K03178/K08844/K03029/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K02138/K03939/K02739/K03937/K02726/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K03942/K02734/K00413/K11351/K02264/K03941/K03063/K03940/K02730/K02728/K02732/K03671
    ## ko05020                                                                                                                                                                                                                                                                    K03968/K03883/K00417/K00416/K03066/K03065/K02127/K00235/K02729/K08738/K03035/K03037/K03030/K03943/K00237/K20858/K03039/K03954/K03033/K00415/K05870/K03953/K07374/K04958/K03934/K02725/K03965/K05863/K03032/K03062/K03064/K03061/K04348/K06268/K09490/K03028/K03097/K10396/K07375/K03880/K03881/K03879/K03083/K04345/K03029/K04962/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K06278/K02138/K03939/K02739/K03937/K02726/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K03942/K02734/K00413/K11351/K02264/K03941/K03063/K03940/K02730/K02728/K02732
    ## ko05010                                                                                                                                                                                                                   K03968/K03883/K00417/K00416/K03066/K03065/K02127/K02183/K00235/K02729/K08738/K03035/K03037/K08957/K03030/K03943/K00237/K02580/K20858/K03039/K04527/K03954/K03033/K00415/K17908/K03953/K07374/K04958/K03934/K02725/K03965/K05863/K03032/K03062/K03064/K03061/K04348/K06268/K08331/K03028/K03097/K08269/K10396/K07375/K03880/K03881/K03879/K03083/K17589/K01408/K03029/K17906/K03031/K11352/K03960/K03950/K06693/K11353/K02134/K03962/K02138/K03939/K02739/K03937/K02726/K08334/K00411/K03963/K03964/K02137/K02727/K03955/K02268/K03966/K00236/K03936/K03036/K02738/K02266/K02263/K03942/K02734/K00413/K11351/K02264/K03941/K03063/K03940/K02730/K02728/K02732/K06704
    ##         Count
    ## ko05014   120
    ## ko05016   100
    ## ko05022   122
    ## ko05012    87
    ## ko05020    85
    ## ko05010    92

And that concludes the tutorial. Feel free to change any of the code to
see how the results and outputs may differ. For example, how do enirched
KO’s differ among upregulated and down regulated genes?
