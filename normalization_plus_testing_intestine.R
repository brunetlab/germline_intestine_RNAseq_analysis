setwd("~/Shuo_Katja")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")
eset <- eset[, pData(eset)$tissue == "intestine"]

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
## 39460  7288
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]

table(fData(eset)$biotype)
##        lincRNA          miRNA          ncRNA          piRNA protein_coding
##             11              1             21              1           7180
##     pseudogene           rRNA         snoRNA          snRNA           tRNA
##             65              6              1              1              1

dim(eset)
## Features  Samples
##     7288       10

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
##  9930412 11289156 12330718 12835614 13284841

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%     50%     75%    100%
## 7256.00 7261.75 7280.00 7285.75 7288.00



dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                                    design = ~ condition + plate,
                                    colData = pData(eset))
dds.deseq <- DESeq(dds)
save(dds.deseq,
     file = "Data/dds_deseq_intestine.RData")

vst <- varianceStabilizingTransformation(dds,
                                         blind = TRUE)
exprs(eset) <- assay(vst)
save(eset,
     file = "Data/eset_vst_intestine.RData")
write.table(exprs(eset),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_intestine.txt")

res <- results(dds.deseq, list("conditionash.2.RNAi", "conditionempty.vector"))
res <- res[order(res$padj),]
sum(res$padj < 0.05, na.rm = TRUE)
## 323
save(res,
     file = "Data/intestine_ash2_vs_control.RData")
write.table(res,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Results/intestine_ash2_vs_control.txt")
