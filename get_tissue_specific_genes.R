#### Instead of testing (tissues are too different to be normalized together), I will follow Wolfgang Huber's advice: "Why not do something pragmatic, like choosing all genes that are higher (lower) in all replicates of gonad than in all replicates of intestine. Or, all genes that are tissue-specific?" and "I would suggest to derive some pragmatic rule for what you call expressed (e.g. RPKM>5) and not expressed (e.g. RPKM<0.5) - you can do various diagnostic plots to tune these parameters."

setwd("~/Shuo_Katja")
library(Biobase)

load("Data/eset_fpkm_intestine_plus_gonad.RData")
eset <- eset.fpkm
## I will only look at empty vector samples:
eset <- eset[, pData(eset)$condition == "empty vector"]


## Get minimum and maximum values per gene and tissue:
min.max <- apply(exprs(eset),
                 1,
                 function(x){
                     min <- tapply(x,
                                   INDEX = pData(eset)$tissue,
                                   FUN = function(y){
                                       return(c(min = min(y),
                                                max = max(y)))
                                   }
                                   )
                 }
                 )



### Difference between minimum of one tissue and maximum of the other tissue must be at least 1 (two-fold) (this way I would get elt-2, for example).

gonad <- names(min.max)[sapply(min.max,
                               function(x){
                                   x$gonad["min"] - x$intestine["max"] >  1
                               })]
length(gonad)
## 5494
is.element(c("pie-1", "mex-5"), gonad)
## TRUE TRUE

intestine <- names(min.max)[sapply(min.max,
                               function(x){
                                   x$intestine["min"] - x$gonad["max"] >  1
                               })]
length(intestine)
## 1418
is.element("elt-2", intestine)
## TRUE


## Write out tables with tissue specific genes and sample-wise FPKM values:

write.table(exprs(eset)[gonad, ],
	col.names = NA,
	sep = "\t",
	quote = FALSE,
	file = "Results/FPKM_gonad_specific_genes.txt")

write.table(exprs(eset)[intestine, ],
        col.names = NA,
        sep = "\t",
        quote =	FALSE,
        file = "Results/FPKM_intestine_specific_genes.txt")
