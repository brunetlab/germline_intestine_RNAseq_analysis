### Although I know that the dispersion shrinking doesn't quite work as good when the data has different degrees of dispersions, I still want to compare gonad vs intestine (empty vector).

setwd("~/Shuo_Katja")
library(Biobase)
library(DESeq2)
library(edgeR)
load("Data/eset_counts.RData")

## Remove the first 4 rows (unmapped reads etc):
eset <- eset[-c(1:4),]

## Filter out genes with low coverage across samples.
## I will only keep genes that have a cpm (counts per million) value of at least 1 in at least 5 samples:

cov.per.sample <- colSums(exprs(eset))
norm.fact <- rep(cov.per.sample, each = nrow(eset))
cpm <- exprs(eset) / norm.fact * 10^6
ind.keep <- rowSums(cpm >= 1) >= 5
table(ind.keep)
## FALSE  TRUE
## 36306 10442
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]


## Instead of CPM, it is better if I go with FPKM values to normalize (to avoid across-sample normalization since tissues are too different):
load("Data/gene_sizes.RData")
gene.sizes <- gene.sizes[fData(eset)$gene.id]
identical(fData(eset)$gene.id, names(gene.sizes))
## TRUE
fData(eset)$gene.size <- gene.sizes

fpkm <- rpkm(exprs(eset),
             gene.length = fData(eset)$gene.size,
             normalized.lib.sizes = FALSE,
             log = FALSE)
eset.fpkm <- eset
exprs(eset.fpkm) <- log2(fpkm + 1)
save(eset.fpkm,
     file = "Data/eset_fpkm_intestine_plus_gonad.RData")

table(fData(eset)$biotype)
##        lincRNA          miRNA          ncRNA          piRNA protein_coding
##             19              1             48              1          10145
##     pseudogene           rRNA         snoRNA          snRNA           tRNA
##            197              6             12             12              1

dim(eset)
## Features  Samples
##    10442       20

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
##  9948852 12578114 13351878 14961878 21954851

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%     50%     75%    100%
## 9159.00  9396.75  9831.50 10156.25 10259.00
