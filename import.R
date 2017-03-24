## Make expressionSet with count matrix.
## module load destiny

setwd("~/Shuo_Katja")
source("~/My_R_functions/make.transparent.R")
library(Biobase)

files <- list.files("Data/Star/Star_output",
                    pattern = "ReadsPerGene",
                    full.names = TRUE)
samples <- gsub("Data/Star/Star_output/", "", files)
samples <- gsub("_ReadsPerGene.out.tab", "", samples)
system(paste("wc -l", files[1]))
## 46752 Data/Star/Star_output/SHPool20-10-CGAGGCTG_S10_ReadsPerGene.out.tab
count.matrix <- matrix(nrow = 46752,
                   ncol = length(files),
                   dimnames = list(NULL, samples))

## counts for unstranded RNA-seq in 2nd column.
for(i in seq(along = files)){
  data.i <- read.table(files[i],
                       colClasses = c("NULL", "integer","NULL", "NULL"),
                       sep = "\t")[,1]
  count.matrix[, samples[i]] <- data.i
  rm(data.i)
}
genes <- read.table(files[i],
                    colClasses = c("character", "NULL","NULL", "NULL"),
                    sep = "\t")[,1]
rownames(count.matrix) <- genes

## Get Sample info:
pData <- read.table("Data/sample_info.txt",
                    sep = "\t",
                    header = TRUE)
pData <- pData[match(colnames(count.matrix), pData$file), ]
pData$group <- paste(pData$tissue, pData$condition)
source("~/My_R_functions/make.transparent.R")
pData$col[pData$tissue == "gonad"
                      & pData$condition == "empty vector"] <- "orangered1"
pData$col[pData$tissue == "gonad"
                      & pData$condition == "ash-2 RNAi"] <- "orangered4"
pData$col[pData$tissue == "intestine"
                      & pData$condition == "empty vector"] <- "deepskyblue3"
pData$col[pData$tissue == "intestine"
                      & pData$condition == "ash-2 RNAi"] <- "midnightblue"
pData$bg <- make.transparent(pData$col, alpha = 200)
pData$dissection <- as.factor(pData$dissection)
pData$plate <- as.factor(pData$plate)
pData$condition <- relevel(as.factor(pData$condition), ref = "empty vector")
pData$tissue <- as.factor(pData$tissue)

## Get gene info:
tab <- read.table("Data/gene_info.txt",
                  sep = ";",
                  stringsAsFactors = FALSE)
fData <- data.frame(gene.id = c(rownames(count.matrix)[1:4], gsub("gene_id ", "", tab$V1)),
                    gene.name = c(rownames(count.matrix)[1:4], gsub(" gene_name ", "", tab$V2)),
                    source = c(rownames(count.matrix)[1:4], gsub(" gene_source ", "", tab$V3)),
                    biotype = c(rownames(count.matrix)[1:4], gsub(" gene_biotype ", "", tab$V4)),
                    stringsAsFactors = FALSE)
rownames(fData) <- fData$gene.name
identical(fData$gene.id, rownames(count.matrix))
## TRUE

eset <- ExpressionSet(assayData = count.matrix)
pData(eset) <- pData
fData(eset) <- fData
rownames(eset) <- fData(eset)$gene.name
colnames(eset) <- paste(pData(eset)$tissue, pData(eset)$condition, pData(eset)$plate)

eset <- eset[, order(colnames(eset))]

save(eset,
     file = "Data/eset_counts.RData")
