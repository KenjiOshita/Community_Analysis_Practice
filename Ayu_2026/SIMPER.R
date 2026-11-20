#-------------------------------------------------------------------------------
# 8.SIMPER（どの属が地点間の違いに効いているか）
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)

Sys.setenv(LANG = "en")

out.dir <- "SIMPER"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 8-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site <- factor(site_names, levels = unique(site_names))

# 8-2. SIMPER（contribution が大きい属ほど、その2地点の違いを作っている）
set.seed(123)
sim <- simper(spe, site, permutations = 99)
summary(sim, ordered = TRUE)

# 8-3. 地点の組み合わせごとに上位5属を保存
top5 <- do.call(rbind, lapply(names(sim), function(nm) {
  d <- head(summary(sim, ordered = TRUE)[[nm]], 5)
  data.frame(Pair = nm, Genus = rownames(d), contribution = d$average, cumsum = d$cumsum)
}))
top5
write.csv(top5, file.path(out.dir, "table_simper_top5.csv"), row.names = FALSE)
