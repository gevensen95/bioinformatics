Phylogenetic Analysis of Doubly Uniparental Inheritance (DUI) in
Bivalves
================
K. Garrett Evensen
5/3/2022

## Introduction to Phylogenetics

We are going to make a phylogenetic tree of bivalves using several
different methods. The first is a maximum parsimony treeusing aligned
cox1 sequences created using MUSCLE. We will create two trees using R
and then two more using RAxML (maximum likelihood) and BEAST (Bayesian).
Once we have learned how to make and graph phylogenetic trees, we can
conduct ancestral state reconstruction of Doubly Uniparental Inheritance
(DUI) in bivalves and determine if DUI is correlated with gonochorism.

Before we can make a phylogenetic tree, we some sequences. There is a
file called cox1.fasta. That file has cytochrome c oxidase sequences for
18 species: 16 bivalves and 2 other mollusks. Next, we need to align
these sequences. I normally use MUSCLE. There are other alignment
softwares that you can use though. You can use MUSCLE by downloading the
software (paste link) and running it on the command line, or use an
online tool (paste links) that uses MUSCLE. If you are using it on the
command line, use the code pasted below:

muscle3.8.31_i86darwin64 -in cox1.fasta -fastaout cox1_align.fasta

Now we can make phylogenetic trees!

If you don’t have treeio or ggtree installed, run this code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("treeio")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.17), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'treeio'

    ## Old packages: 'BiocManager', 'openssl', 'tinytex'

``` r
BiocManager::install("ggtree")
```

    ## Bioconductor version 3.15 (BiocManager 1.30.17), R 4.2.0 (2022-04-22)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'ggtree'

    ## Old packages: 'BiocManager', 'openssl', 'tinytex'

``` r
require(ape)
```

    ## Loading required package: ape

``` r
require(phangorn)
```

    ## Loading required package: phangorn

``` r
require(phytools)
```

    ## Loading required package: phytools

    ## Loading required package: maps

``` r
require(ggtree)
```

    ## Loading required package: ggtree

    ## ggtree v3.4.0 For help: https://yulab-smu.top/treedata-book/
    ## 
    ## If you use the ggtree package suite in published research, please cite
    ## the appropriate paper(s):
    ## 
    ## Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
    ## ggtree: an R package for visualization and annotation of phylogenetic
    ## trees with their covariates and other associated data. Methods in
    ## Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
    ## 
    ## LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
    ## Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
    ## for phylogenetic tree input and output with richly annotated and
    ## associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
    ## doi: 10.1093/molbev/msz240
    ## 
    ## Guangchuang Yu. Using ggtree to visualize data on tree-like structures.
    ## Current Protocols in Bioinformatics. 2020, 69:e96. doi:10.1002/cpbi.96

    ## 
    ## Attaching package: 'ggtree'

    ## The following object is masked from 'package:ape':
    ## 
    ##     rotate

``` r
require(treeio)
```

    ## Loading required package: treeio

    ## treeio v1.20.0 For help: https://yulab-smu.top/treedata-book/
    ## 
    ## If you use the ggtree package suite in published research, please cite
    ## the appropriate paper(s):
    ## 
    ## LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
    ## Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
    ## for phylogenetic tree input and output with richly annotated and
    ## associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
    ## doi: 10.1093/molbev/msz240
    ## 
    ## Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
    ## ggtree: an R package for visualization and annotation of phylogenetic
    ## trees with their covariates and other associated data. Methods in
    ## Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
    ## 
    ## S Xu, Z Dai, P Guo, X Fu, S Liu, L Zhou, W Tang, T Feng, M Chen, L
    ## Zhan, T Wu, E Hu, Y Jiang, X Bo, G Yu. ggtreeExtra: Compact
    ## visualization of richly annotated phylogenetic data. Molecular Biology
    ## and Evolution. 2021, 38(9):4039-4042. doi: 10.1093/molbev/msab166

    ## 
    ## Attaching package: 'treeio'

    ## The following object is masked from 'package:phytools':
    ## 
    ##     read.newick

    ## The following object is masked from 'package:ape':
    ## 
    ##     drop.tip

First we will import the aligned sequences in fasta format and then
convert it to a phyDat object.

``` r
cox1.dna <- read.dna("cox1_align.fasta", format="fasta")
cox1.phyDat <- as.phyDat(cox1.dna)
```

# Maximum Parsimony

The first tree we will make is a maximum parsimony (MP) tree (think
Occam’s Razor).

``` r
mp.tree<-pratchet(cox1.phyDat) 
```

    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"
    ## [1] "Best pscore so far: 5783"

``` r
mp.tree$edge.length<-runif(n=nrow(mp.tree$edge))
plot.phylo(mp.tree)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Maximum Likelihood

The next method is maximum likelihood (ML). Before we can make the tree,
we have to determine which is the best model for making an ML tree.

``` r
mt.result<-modelTest(pml(mp.tree,data=cox1.phyDat)) #determine which is the best model for ML
```

    ## [1] "JC+I"
    ## [1] "JC+G"
    ## [1] "JC+G+I"
    ## [1] "F81+I"
    ## [1] "F81+G"
    ## [1] "F81+G+I"
    ## [1] "K80+I"
    ## [1] "K80+G"
    ## [1] "K80+G+I"
    ## [1] "HKY+I"
    ## [1] "HKY+G"
    ## [1] "HKY+G+I"
    ## [1] "SYM+I"
    ## [1] "SYM+G"
    ## [1] "SYM+G+I"
    ## [1] "GTR+I"
    ## [1] "GTR+G"
    ## [1] "GTR+G+I"

``` r
rnd<-function(x,...) if(is.numeric(x)) round(x,...) else x
as.data.frame(lapply(mt.result,rnd,digits=2))
```

    ##      Model df    logLik      AIC AICw     AICc AICcw      BIC
    ## 1       JC 45 -25225.49 50540.99    0 50543.32     0 50788.83
    ## 2     JC+I 46 -23871.99 47835.99    0 47838.43     0 48089.34
    ## 3     JC+G 46 -23503.27 47098.54    0 47100.98     0 47351.89
    ## 4   JC+G+I 47 -23447.19 46988.39    0 46990.93     0 47247.25
    ## 5      F81 48 -25069.80 50235.60    0 50238.25     0 50499.97
    ## 6    F81+I 49 -23679.66 47457.32    0 47460.09     0 47727.20
    ## 7    F81+G 49 -23146.45 46390.91    0 46393.67     0 46660.78
    ## 8  F81+G+I 50 -23126.08 46352.17    0 46355.05     0 46627.55
    ## 9      K80 46 -25016.08 50124.15    0 50126.59     0 50377.51
    ## 10   K80+I 47 -23642.95 47379.89    0 47382.43     0 47638.75
    ## 11   K80+G 47 -23193.26 46480.52    0 46483.06     0 46739.38
    ## 12 K80+G+I 48 -23152.52 46401.05    0 46403.70     0 46665.42
    ## 13     HKY 49 -24755.34 49608.68    0 49611.45     0 49878.56
    ## 14   HKY+I 50 -23315.51 46731.02    0 46733.90     0 47006.40
    ## 15   HKY+G 50 -22586.31 45272.63    0 45275.51     0 45548.01
    ## 16 HKY+G+I 51 -22568.50 45239.00    0 45241.99     0 45519.89
    ## 17     SYM 50 -24730.11 49560.21    0 49563.09     0 49835.60
    ## 18   SYM+I 51 -23374.22 46850.43    0 46853.43     0 47131.33
    ## 19   SYM+G 51 -23021.58 46145.16    0 46148.16     0 46426.06
    ## 20 SYM+G+I 52 -22957.16 46018.32    0 46021.44     0 46304.72
    ## 21     GTR 53 -24612.93 49331.86    0 49335.10     0 49623.77
    ## 22   GTR+I 54 -23222.67 46553.35    0 46556.71     0 46850.76
    ## 23   GTR+G 54 -22497.68 45103.36    0 45106.72     0 45400.78
    ## 24 GTR+G+I 55 -22486.79 45083.58    1 45087.07     1 45386.50

Based on these results, the GTR+G+I model is the best model for ML. In
the optim.pml function, we set optGamma=TRUE and optInv=TRUE to add the
G and I to the GTR model.

``` r
obj<-pml(tree=rtree(n=length(cox1.phyDat),tip.label=names(cox1.phyDat),
                    rooted=FALSE),data=cox1.phyDat)
fit.ml <- optim.pml(obj,optEdge=TRUE,model="GTR",optGamma=TRUE, optInv=TRUE,rearrangement="NNI")
```

    ## only one rate class, ignored optGamma

    ## optimize edge weights:  -35225.61 --> -27920.98 
    ## optimize base frequencies:  -27920.98 --> -27726.94 
    ## optimize rate matrix:  -27726.94 --> -27375.93 
    ## optimize invariant sites:  -27375.93 --> -25337.28 
    ## optimize edge weights:  -25337.28 --> -25272.62 
    ##  optimize topology:  -25272.62 --> -25171.94 
    ##  optimize topology:  -25171.94 --> -25122.08 
    ##  optimize topology:  -25122.08 --> -25073.49 
    ##  optimize topology:  -25073.49 --> -24989.34 
    ##  optimize topology:  -24989.34 --> -24964.42 
    ## NNI moves:  16 
    ## optimize base frequencies:  -24964.42 --> -24954.3 
    ## optimize rate matrix:  -24954.3 --> -24925.58 
    ## optimize invariant sites:  -24925.58 --> -24920.02 
    ## optimize edge weights:  -24920.02 --> -24908.04 
    ##  optimize topology:  -24908.04 --> -24852.95 
    ##  optimize topology:  -24852.95 --> -24790.56 
    ##  optimize topology:  -24790.56 --> -24641.64 
    ##  optimize topology:  -24641.64 --> -24601.19 
    ##  optimize topology:  -24601.19 --> -24567.09 
    ## NNI moves:  11 
    ## optimize base frequencies:  -24567.09 --> -24563.84 
    ## optimize rate matrix:  -24563.84 --> -24561.93 
    ## optimize invariant sites:  -24561.93 --> -24555.45 
    ## optimize edge weights:  -24555.45 --> -24550.91 
    ##  optimize topology:  -24550.91 --> -24488.37 
    ##  optimize topology:  -24488.37 --> -24310.01 
    ##  optimize topology:  -24310.01 --> -24158.99 
    ##  optimize topology:  -24158.99 --> -24150.8 
    ##  optimize topology:  -24150.8 --> -24136.36 
    ## NNI moves:  12 
    ## optimize base frequencies:  -24136.36 --> -24132.33 
    ## optimize rate matrix:  -24132.33 --> -24130.01 
    ## optimize invariant sites:  -24130.01 --> -24126.24 
    ## optimize edge weights:  -24126.24 --> -24124.72 
    ##  optimize topology:  -24124.72 --> -23974.71 
    ##  optimize topology:  -23974.71 --> -23899.42 
    ##  optimize topology:  -23899.42 --> -23778.2 
    ##  optimize topology:  -23778.2 --> -23714.7 
    ##  optimize topology:  -23714.7 --> -23619.64 
    ## NNI moves:  6 
    ## optimize base frequencies:  -23619.64 --> -23617.67 
    ## optimize rate matrix:  -23617.67 --> -23615.78 
    ## optimize invariant sites:  -23615.78 --> -23614.7 
    ## optimize edge weights:  -23614.7 --> -23613.99 
    ##  optimize topology:  -23613.99 --> -23511.63 
    ##  optimize topology:  -23511.63 --> -23416.07 
    ##  optimize topology:  -23416.07 --> -23386.17 
    ##  optimize topology:  -23386.17 --> -23308.52 
    ##  optimize topology:  -23308.52 --> -23288.58 
    ## NNI moves:  8 
    ## optimize base frequencies:  -23288.58 --> -23287.49 
    ## optimize rate matrix:  -23287.49 --> -23285.9 
    ## optimize invariant sites:  -23285.9 --> -23285.76 
    ## optimize edge weights:  -23285.76 --> -23285.56 
    ##  optimize topology:  -23285.56 --> -23266.99 
    ##  optimize topology:  -23266.99 --> -23257.29 
    ##  optimize topology:  -23257.29 --> -23250.92 
    ##  optimize topology:  -23250.92 --> -23232.79 
    ##  optimize topology:  -23232.79 --> -23223.41 
    ## NNI moves:  6 
    ## optimize base frequencies:  -23223.41 --> -23223.31 
    ## optimize rate matrix:  -23223.31 --> -23223.03 
    ## optimize invariant sites:  -23223.03 --> -23222.99 
    ## optimize edge weights:  -23222.99 --> -23222.93 
    ##  optimize topology:  -23222.93 --> -23222.05 
    ##  optimize topology:  -23222.05 --> -23216.91 
    ##  optimize topology:  -23216.91 --> -23216.91 
    ## NNI moves:  2 
    ## optimize base frequencies:  -23216.91 --> -23216.88 
    ## optimize rate matrix:  -23216.88 --> -23216.87 
    ## optimize invariant sites:  -23216.87 --> -23216.86 
    ## optimize edge weights:  -23216.86 --> -23216.84 
    ##  optimize topology:  -23216.84 --> -23216.84 
    ## NNI moves:  0 
    ## optimize base frequencies:  -23216.84 --> -23216.84 
    ## optimize rate matrix:  -23216.84 --> -23216.84 
    ## optimize invariant sites:  -23216.84 --> -23216.84 
    ## optimize edge weights:  -23216.84 --> -23216.83 
    ## optimize base frequencies:  -23216.83 --> -23216.83 
    ## optimize rate matrix:  -23216.83 --> -23216.83 
    ## optimize invariant sites:  -23216.83 --> -23216.83 
    ## optimize edge weights:  -23216.83 --> -23216.83 
    ## optimize base frequencies:  -23216.83 --> -23216.83 
    ## optimize rate matrix:  -23216.83 --> -23216.83 
    ## optimize invariant sites:  -23216.83 --> -23216.83 
    ## optimize edge weights:  -23216.83 --> -23216.83 
    ## optimize base frequencies:  -23216.83 --> -23216.83 
    ## optimize rate matrix:  -23216.83 --> -23216.83 
    ## optimize invariant sites:  -23216.83 --> -23216.83 
    ## optimize edge weights:  -23216.83 --> -23216.83

``` r
plot.phylo(midpoint.root(fit.ml$tree))
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

How does the ML tree compare to the MP tree? Are they the same?

# Bootstrapping Maximum Likelihood

To improve our confidence in the strenght of nodes, we will use
bootstrapping. You can use 100, but 1000 is better. Here, we will only
use 100 for the sake of computational time. For plotting, we will keep
only bootstrap values greater than 50. A value of 70 is roughly
equivalent to a p value of 0.05.

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Maximum Likelihood with RAxML

To make a tree using RAxML
(<https://cme.h-its.org/exelixis/web/software/raxml/>), you can use the
code provided on Github, or use the code copied below:

raxmlHPC-PTHREADS-AVX -s cox1_align.fasta -f a -k -# 1000 -n cox1 -m
GTRGAMMAI -p 912046 -x 14826 -T 2

This code bootstraps 1000 times and uses the GTR+G+I model. Make sure to
read the vignette to understand all the possible options for RAxML.
There will be five different output files. We will use the one with
bipartitionsBranchLabels. You can draw the phylogenetic tree in two
different ways. The first using plot.phylo in the ape package, and the
second uses ggtree.

``` r
cox1_raxml <- read.raxml('RAxML_bipartitionsBranchLabels.cox1')
```

``` r
plot.phylo(cox1_raxml@phylo)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Next, we will use ggtree. This package gives us lots of different
options for graphing. It works similar to ggplot. For now, we will keep
it simple.

``` r
ggtree(cox1_raxml) + xlim(0, 4) +
  geom_tiplab(align = TRUE) + 
  geom_text2(aes(label=bootstrap), hjust = 1, vjust = -0.1)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

BEAST

Next, we are going to use BEAST2 (<https://www.beast2.org>) to make a
Bayesian tree using the same alignment file. Make sure to read the full
instructions online for proper installation. Below I’ve pasted quick
instructions on what I did to make this tree.

Convert alignment to nexus format:
<http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi>

Open BEAUti, load the .nex file by clicking the plus on the bottom left
of the window. Under Site Model set the substitution model to GTR. It
should be JC69 by default. Then generate the BEAST .xml file by clicking
file and saving the .xml file. It will defualt to 10,000,000
generations, but you can change this if you want under MCMC.

Next, open BEAST choose the .xml you just generated. I called mine
cox1.xml because we are working with the cox1 gene. You can keep the
Thread Pool Size set to automatic, or you can set the number of threads
depending on your machine and if you are doing anything else while BEAST
is running. Check Use BEAGLE library if available. Then click Run.

Finally, open TreeAnnotator and load the .trees file. Set the burnin to
be 25%. Name the output and click Run. I’ve named the file
cox1_beast.tre. Now we can import the BEAST tree and graph it.

``` r
cox1_beast <- read.beast('cox1_beast.tre')
```

``` r
ggtree(cox1_beast) + xlim(0, 3) +
  geom_tiplab(align = TRUE) + 
  geom_text2(aes(label=round(as.numeric(posterior), 2), x=branch), hjust = 0.75, vjust=0)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# Compare RAxML and BEAST Trees

Frequently, you will want to compare the topologies of trees made using
the same alignment file, but with two different methods. We can do this
using the phytools package.

``` r
t1 <- force.ultrametric(cox1_raxml@phylo, "nnls")
```

    ## ***************************************************************
    ## *                          Note:                              *
    ## *    force.ultrametric does not include a formal method to    *
    ## *    ultrametricize a tree & should only be used to coerce    *
    ## *   a phylogeny that fails is.ultramtric due to rounding --   *
    ## *    not as a substitute for formal rate-smoothing methods.   *
    ## ***************************************************************

``` r
t2 <- force.ultrametric(cox1_beast@phylo, "nnls")
```

    ## ***************************************************************
    ## *                          Note:                              *
    ## *    force.ultrametric does not include a formal method to    *
    ## *    ultrametricize a tree & should only be used to coerce    *
    ## *   a phylogeny that fails is.ultramtric due to rounding --   *
    ## *    not as a substitute for formal rate-smoothing methods.   *
    ## ***************************************************************

``` r
obj2 <- cophylo(t1, t2)
```

    ## Rotating nodes to optimize matching...
    ## Done.

``` r
plot(obj2,link.type="curved",link.lwd=4,link.lty="solid",
     link.col=make.transparent("blue",0.4),fsize=1)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

As you can see, the topologies of the two trees are nearly the same. Are
they the same as the MP tree and first ML tree we made? How are they
different?

# Rooting a BEAST Tree

The best way to root a tree with BEAST is to do it when you are first
generating the BEAST file. Remake the BEAST file, but let’s add an
outgroup this time. Set up the file as we did the first time. Then click
Priors, then click Add prior, and move all the taxa to the right side,
except Conus loroisii and Nautilus pompilius. These taxa will from a
monophyletic outgroup. C. loroisii and N. pompilius are members of
Gastropoda and Cephalopoda, respectively. Click OK and save the file as
normal. Run BEAST as we did earler.

Let’s import the new tree and graph it. We will add extra information on
the tree, such as if a species is gonochoristic/hermaphroditic, exhibits
DUI/No DUI, or is from a particular clade. To identify node numbers, add
geom_text(aes(label=node))

``` r
cox1_beast_root <- read.beast('cox1_beast_root.tre')

c.names <- read.csv(file = "names.csv", header=FALSE) #import an csv file of the original tip labs and what we want to rename them as
cox1_beast_root <- rename_taxa(cox1_beast_root, c.names) #rename the taxa
bvData<-read.csv("bvData.csv")

ggtree(cox1_beast_root) %<+% bvData + xlim(0, 4) +
  geom_tiplab(align = TRUE) + 
  geom_cladelab(node=38, label="Heterodonta", align=TRUE,  
                offset = 1.7) +
  geom_cladelab(node=43, label="Paleoheterodonta", align=TRUE,  
                offset = 1.7) +
  geom_cladelab(node=32, label="Pteriomorphia", align=TRUE,  
                offset = 1.7) +
  geom_cladelab(node=26, label="Outgroup", align=TRUE,  
                offset = 1.7) +
  geom_tippoint(aes(color = DUI, shape = mode)) +
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), hjust = 0.75, vjust=0) 
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

That covers the basics of making and graphing trees. Congrats! To test
yourself, change any of the inputs in the ggtree object, or even try
making a tree of your own using different sequences, such as 16S rRNA or
even a gene tree.

## Ancestral State Reconstruction and Correlation of Discrete Characters

An important step of understanding a discrete character in a taxonomic
group is determining when the trait evolved. Furthermore, many discrete
characters are correlated with other, possibly related traits. In this
particular case in bivalves, we want to know if gonochorism and no DUI
are ancestral states of bivalves.

# Ancestral State Reconstruction

We will use ancestral state reconstruction of
gonochorism/hermaphroditism and DUI/No DUI using the phytools package by
Liam Revell. We will do this using the rooted BEAST tree we made
earlier. First we will evaluate if gonochorism or hermaphroditism is the
ancestral state of bivalves.

``` r
sex_modes <- as.factor(setNames(bvData$mode, bvData$Bivalve))

bivalve.trees<-make.simmap(drop.tip(cox1_beast_root@phylo,setdiff(cox1_beast_root@phylo$tip.label,names(sex_modes))),sex_modes,model="SYM",nsim=100)
```

    ## make.simmap is sampling character histories conditioned on
    ## the transition matrix
    ## 
    ## Q =
    ##                Gonochoristic Hermaphroditic
    ## Gonochoristic      -1.400628       1.400628
    ## Hermaphroditic      1.400628      -1.400628
    ## (estimated using likelihood);
    ## and (mean) root node prior probabilities
    ## pi =
    ##  Gonochoristic Hermaphroditic 
    ##            0.5            0.5

    ## Done.

``` r
obj <- summary(bivalve.trees, plot = FALSE)
cols<-setNames(palette()[1:2],levels(sex_modes))
plot(obj,colors=cols,fsize=0.8,cex=c(0.5,0.3))
add.simmap.legend(colors=cols, prompt = FALSE, y=2)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Now, we will do the same for DUI.

``` r
DUI_modes <- as.factor(setNames(bvData$DUI, bvData$Bivalve))

bivalve.trees<-make.simmap(drop.tip(cox1_beast_root@phylo,setdiff(cox1_beast_root@phylo$tip.label,names(DUI_modes))),DUI_modes,model="SYM",nsim=100)
```

    ## make.simmap is sampling character histories conditioned on
    ## the transition matrix
    ## 
    ## Q =
    ##              DUI    No DUI
    ## DUI    -38.88773  38.88773
    ## No DUI  38.88773 -38.88773
    ## (estimated using likelihood);
    ## and (mean) root node prior probabilities
    ## pi =
    ##    DUI No DUI 
    ##    0.5    0.5

    ## Done.

``` r
obj <- summary(bivalve.trees, plot = FALSE)
cols<-setNames(palette()[1:2],levels(DUI_modes))
plot(obj,colors=cols,fsize=0.8,cex=c(0.5,0.3))
add.simmap.legend(colors=cols, prompt = FALSE, y = 2)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Based on these two figures, we can conclude that the ancestral bivalve
was gonochoristic and likely did not exhibit DUI. Based on stochastic
mapping of DUI, this phenomenon likely evolved three separate times in
bivalves, once in each bivalve subclass.

# Correlation of Discrete Traits

``` r
par(mfrow=c(1,2))
plot(force.ultrametric(cox1_beast_root@phylo),show.tip.label=FALSE,no.margin=TRUE)
```

    ## ***************************************************************
    ## *                          Note:                              *
    ## *    force.ultrametric does not include a formal method to    *
    ## *    ultrametricize a tree & should only be used to coerce    *
    ## *   a phylogeny that fails is.ultramtric due to rounding --   *
    ## *    not as a substitute for formal rate-smoothing methods.   *
    ## ***************************************************************

``` r
par(fg="transparent")
tiplabels(pie=to.matrix(sex_modes[cox1_beast_root@phylo$tip.label],c("Gonochoristic","Hermaphroditic")),piecol=c("red","blue"),cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("red","blue"),c("Gonochoristic","Hermaphroditic")),prompt=FALSE, x=0.1,y=7,fsize=1)
plot(force.ultrametric(cox1_beast_root@phylo),show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
```

    ## ***************************************************************
    ## *                          Note:                              *
    ## *    force.ultrametric does not include a formal method to    *
    ## *    ultrametricize a tree & should only be used to coerce    *
    ## *   a phylogeny that fails is.ultramtric due to rounding --   *
    ## *    not as a substitute for formal rate-smoothing methods.   *
    ## ***************************************************************

``` r
tiplabels(pie=to.matrix(DUI_modes[cox1_beast_root@phylo$tip.label],c("DUI","No DUI")),piecol=c("brown","cyan"),cex=0.3)
add.simmap.legend(colors=setNames(c("brown","cyan"),c("DUI","No DUI")),prompt=FALSE, x=0.7,y=7,fsize=1)
```

![](phylogenetics_intro_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Here, we want to match up the colors. Are red circles (gonochorism)
matched up with brown circles (DUI)? Are blue circles (hermaphroditism)
matched up with cyan circles (No DUI)?

Now, we use a model to determine if sex mode depends on DUI mode and
vice versa.

``` r
fit.sd <- fitPagel(cox1_beast_root@phylo, sex_modes, DUI_modes)
fit.sd
```

    ## 
    ## Pagel's binary character correlation test:
    ## 
    ## Assumes "ARD" substitution model for both characters
    ## 
    ## Independent model rate matrix:
    ##                       Gonochoristic|DUI Gonochoristic|No DUI Hermaphroditic|DUI
    ## Gonochoristic|DUI            -15.350659             6.229934           9.120725
    ## Gonochoristic|No DUI           4.178212           -13.298937           0.000000
    ## Hermaphroditic|DUI            27.344458             0.000000         -33.574392
    ## Hermaphroditic|No DUI          0.000000            27.344458           4.178212
    ##                       Hermaphroditic|No DUI
    ## Gonochoristic|DUI                  0.000000
    ## Gonochoristic|No DUI               9.120725
    ## Hermaphroditic|DUI                 6.229934
    ## Hermaphroditic|No DUI            -31.522670
    ## 
    ## Dependent (x & y) model rate matrix:
    ##                       Gonochoristic|DUI Gonochoristic|No DUI Hermaphroditic|DUI
    ## Gonochoristic|DUI             -38.88773             38.88773           0.000000
    ## Gonochoristic|No DUI           38.77820            -49.91305           0.000000
    ## Hermaphroditic|DUI              0.00000              0.00000         -11.162863
    ## Hermaphroditic|No DUI           0.00000             19.85284           2.337117
    ##                       Hermaphroditic|No DUI
    ## Gonochoristic|DUI                   0.00000
    ## Gonochoristic|No DUI               11.13485
    ## Hermaphroditic|DUI                 11.16286
    ## Hermaphroditic|No DUI             -22.18996
    ## 
    ## Model fit:
    ##             log-likelihood      AIC
    ## independent      -29.78635 67.57271
    ## dependent        -28.61224 73.22447
    ## 
    ## Hypothesis test result:
    ##   likelihood-ratio:  2.34823 
    ##   p-value:  0.672002 
    ## 
    ## Model fitting method used was fitMk

Based on these results, we cannot conclude that sex_mode and DUI_mode
are correlated because the p value is greater than 0.05. However, there
is some evidence to support that gonochorism and DUI are correlated as
the p value is approximately 0.1. How would the results change if the
RAxML tree was used, or if a completely new/larger tree was used? It is
important to note that every bivalve species that exhibits DUI is a
gonochoristic species, except Semimytilus algosus, which is
hermaphroditic.
