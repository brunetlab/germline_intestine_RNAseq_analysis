setwd("~/Shuo_Katja")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")
eset <- eset[, pData(eset)$tissue == "gonad"]

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
## 37675  9073
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]

table(fData(eset)$biotype)
##        lincRNA          ncRNA protein_coding     pseudogene           rRNA
##             11             37           8833            162              6
##         snoRNA          snRNA           tRNA
##             12             11              1

dim(eset)
## Features  Samples
##     9073       10

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
## 13402501 14437145 14978755 15969169 21941589

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%     50%     75%    100%
## 8981.00 9053.75 9069.50 9072.75 9073.00



dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                                    design = ~ condition + plate,
                                    colData = pData(eset))
dds.deseq <- DESeq(dds)
save(dds.deseq,
     file = "Data/dds_deseq_gonad.RData")

vst <- varianceStabilizingTransformation(dds,
                                         blind = TRUE)
exprs(eset) <- assay(vst)
save(eset,
     file = "Data/eset_vst_gonad.RData")
write.table(exprs(eset),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_gonad.txt")

res <- results(dds.deseq, list("conditionash.2.RNAi", "conditionempty.vector"))
res <- res[order(res$padj),]
sum(res$padj < 0.05, na.rm = TRUE)
## 12 (fold changes are all very small, top 20: nothing reaches 0.5)
save(res,
     file = "Data/gonad_ash2_vs_control.RData")
write.table(res,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Results/gonad_ash2_vs_control.txt")
