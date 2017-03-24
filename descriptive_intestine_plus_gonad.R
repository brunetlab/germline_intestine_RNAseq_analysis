setwd("~/Shuo_Katja")
library(Biobase)
library(RColorBrewer)

load("Data/eset_fpkm_intestine_plus_gonad.RData")
eset <- eset.fpkm

pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:4]
legend <- unique(pData(eset)[, c("group", "col", "bg")])

pdf("Results/pca_fpkm_intestine_plus_gonad.pdf",
    width = 5, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 10.1),
    xpd = TRUE)
plot(pca$x[, "PC1"],
     pca$x[, "PC2"],
     pch = 21,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""))
legend(120, 60,
       col = legend$col,
       pt.bg = legend$bg,
       pch = 21,
       pt.cex = 1.3,
       legend = legend$group,
       bty = "n")
dev.off()



### It looks like intestine doesn't express as many genes and/or gonad expresses some genes very highly. What is it?

## gene-wise median over all samples within one tissue:

avg.intestine <- apply(exprs(eset)[, pData(eset)$tissue == "intestine"],
                       1,
                       median)
avg.gonad <- apply(exprs(eset)[, pData(eset)$tissue == "gonad"],
                   1,
                   median)

sum(avg.gonad > avg.intestine) / length(avg.gonad)
## --> 0.763% of genes higher expressed in gonad

sum(avg.gonad == 0)
## 185 in gonad with zero read counts
sum(avg.intestine == 0)
## 682 in intestine with zero read counts

quantile(avg.gonad, c(0.9, 0.95, 0.99))
##      90%      95%      99%
## 6.580685 7.477621 9.774601

quantile(avg.intestine, c(0.9, 0.95, 0.99))
##      90%      95%      99%
## 4.914382 6.338857 9.031434

