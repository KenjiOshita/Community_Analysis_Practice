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

# 3-2. Chord distance（本番と同じ：行を正規化 -> ユークリッド距離）
spe.norm <- decostand(spe, "normalize")
chord    <- vegdist(spe.norm, method = "euclidean")

# 3-3. 手法の比較（cophenetic correlation：1に近いほど元の距離を再現）
methods   <- c("single", "complete", "average", "ward.D2")   # average = UPGMA
coph_cor  <- sapply(methods, function(m) cor(chord, cophenetic(hclust(chord, method = m))))
round(coph_cor, 3)


# 3-4. クラスター数の決定
library(cluster)

hc <- hclust(chord, method = "average")   # UPGMA
n  <- nrow(spe)

# (a) 平均シルエット幅（大きいほど良い）
Si <- numeric(n)
for (k.i in 2:(n - 1)) {
  sil.i <- silhouette(cutree(hc, k = k.i), chord)
  Si[k.i] <- summary(sil.i)$avg.width
}
k.best.sil <- which.max(Si)

# (b) 行列相関（Pearson）：クラスター所属から作った距離と元の距離の相関
grpdist <- function(X) {
  gr <- as.data.frame(as.factor(X))
  daisy(gr, "gower")
}
kt <- data.frame(k = 1:n, r = 0)
for (i in 2:(n - 1)) {
  kt[i, 2] <- cor(chord, grpdist(cutree(hc, i)), method = "pearson")
}
k.best.cor <- which.max(kt$r)

# 作図
draw_optimal_k <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 4))
  plot(1:n, Si, type = "h",
       main = "Silhouette-optimal number of clusters",
       xlab = "k (number of clusters)", ylab = "Average silhouette width")
  axis(1, k.best.sil, paste("optimum", k.best.sil, sep = "\n"),
       col = "red", font = 2, col.axis = "red")
  points(k.best.sil, max(Si), pch = 16, col = "red", cex = 1.5)
  
  plot(kt$k, kt$r, type = "h",
       main = "Matrix correlation-optimal number of clusters",
       xlab = "k (number of clusters)", ylab = "Pearson's correlation")
  axis(1, k.best.cor, paste("optimum", k.best.cor, sep = "\n"),
       col = "red", font = 2, col.axis = "red")
  points(k.best.cor, max(kt$r), pch = 16, col = "red", cex = 1.5)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}
draw_optimal_k()

png(file.path(out.dir, "fig3_optimal_k.png"), width = 2600, height = 1300, res = 300)
draw_optimal_k()
dev.off()

write.csv(data.frame(k = kt$k, avg_silhouette_width = Si, pearson_r = kt$r),
          file.path(out.dir, "table1_optimal_k.csv"), row.names = FALSE)

# 使うクラスター数（本番と同じく、両方の結果を見て自分で決める）
k.best.sil; k.best.cor
k <- k.best.sil   # <- 必要なら手動で変更（例: k <- 4）


# 3-5. UPGMA dendrogram（ラベルを地点の色で塗る + クラスターの枠）
dend <- as.dendrogram(hc)
labels_colors(dend) <- unname(site_cols[sub("_[0-9]+$", "", labels(dend))])
dend <- set(dend, "labels_cex", 0.8)
dend <- set(dend, "branches_k_color", k = k)   # 決めた k で枝を色分け

draw_dend <- function() {
  par(mar = c(5, 4, 3, 1))
  plot(dend, ylab = "Chord distance",
       main = paste0("UPGMA, k = ", k,
                     " (cophenetic r = ", round(coph_cor["average"], 3), ")"))
  legend("topright", legend = site_order, col = site_cols, pch = 15, bty = "n")
}
draw_dend()

png(file.path(out.dir, "fig4_dendrogram_UPGMA.png"), width = 2400, height = 1600, res = 300)
draw_dend()
dev.off()