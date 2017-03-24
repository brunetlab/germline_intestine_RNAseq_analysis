setwd("~/Shuo_Katja")
library(Biobase)
library(RColorBrewer)

load("Data/eset_vst_intestine.RData")

pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:4]
legend <- unique(pData(eset)[, c("group", "col", "bg")])

pdf("Results/pca_vst_intestine.pdf",
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
     xlim = c(-80, 60),
     ylim = c(-60, 60),
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""))
legend(65, 60,
       col = legend$col,
       pt.bg = legend$bg,
       pch = 21,
       pt.cex = 1.3,
       legend = legend$group,
       bty = "n")
dev.off()



