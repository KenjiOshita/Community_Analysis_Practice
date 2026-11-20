#-------------------------------------------------------------------------------
# 3.Clustering (UPGMA)
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(dendextend)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "Clustering"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 3-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)   # 地点の順番を変えたいとき: c("Upper","Middle",...)
site_cols  <- setNames(brewer.pal(8, "Set2")[(seq_along(site_order) - 1) %% 8 + 1], site_order)

# 3-2. Bray-Curtis dissimilarity
bray <- vegdist(spe, method = "bray")

# 3-3. 手法の比較（cophenetic correlation：1に近いほど元の距離を再現）
methods   <- c("single", "complete", "average", "ward.D2")   # average = UPGMA
coph_cor  <- sapply(methods, function(m) cor(bray, cophenetic(hclust(bray, method = m))))
round(coph_cor, 3)

# 3-4. UPGMA dendrogram（ラベルを地点の色で塗る）
hc   <- hclust(bray, method = "average")
dend <- as.dendrogram(hc)
labels_colors(dend) <- unname(site_cols[sub("_[0-9]+$", "", labels(dend))])
dend <- set(dend, "labels_cex", 0.8)

plot(dend, ylab = "Bray-Curtis dissimilarity",
     main = paste("UPGMA (cophenetic r =", round(coph_cor["average"], 3), ")"))
legend("topright", legend = site_order, col = site_cols, pch = 15, bty = "n")

png(file.path(out.dir, "fig3_dendrogram_UPGMA.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 4, 3, 1))
plot(dend, ylab = "Bray-Curtis dissimilarity",
     main = paste("UPGMA (cophenetic r =", round(coph_cor["average"], 3), ")"))
legend("topright", legend = site_order, col = site_cols, pch = 15, bty = "n")
dev.off()

