#-------------------------------------------------------------------------------
# 7.Pairwise PERMANOVA（どの地点とどの地点が違うか）
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "PERMANOVA"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 7-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)
site <- factor(site_names, levels = site_order)

bray <- vegdist(spe, method = "bray")

# 7-2. 2地点ずつ取り出して PERMANOVA
set.seed(123)
pairs <- combn(site_order, 2)
pw <- data.frame()
for (i in 1:ncol(pairs)) {
  keep <- site %in% pairs[, i]
  d    <- as.dist(as.matrix(bray)[keep, keep])
  grp  <- droplevels(site[keep])
  res  <- adonis2(d ~ grp, permutations = 999)
  pw   <- rbind(pw, data.frame(Pair = paste(pairs[, i], collapse = " vs "),
                               R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1]))
}
pw$p_adj <- p.adjust(pw$p, method = "holm")   # 比較回数が多いので補正
pw
write.csv(pw, file.path(out.dir, "table_pairwise_permanova.csv"), row.names = FALSE)

