#-------------------------------------------------------------------------------
# 5.PERMANOVA and PERMDISP
#-------------------------------------------------------------------------------
# PERMANOVA: 地点ごとの「平均的な組成」に差があるか
# PERMDISP : 地点ごとの「ばらつき」に差があるか（PERMANOVAとセットで見る）
rm(list = ls())

library(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "PERMANOVA"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 5-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)
site <- factor(site_names, levels = site_order)
site_cols <- setNames(brewer.pal(8, "Set2")[(seq_along(site_order) - 1) %% 8 + 1], site_order)

bray <- vegdist(spe, method = "bray")
jac  <- vegdist(spe, method = "jaccard", binary = TRUE)

# 5-2. PERMANOVA（R2 = 地点で説明できる割合、Pr(>F) = p値）
set.seed(123)   # 並べ替え検定は乱数を使うので固定
perm_bray <- adonis2(bray ~ site, permutations = 999)
perm_jac  <- adonis2(jac  ~ site, permutations = 999)
perm_bray
perm_jac
write.csv(as.data.frame(perm_bray), file.path(out.dir, "table_permanova_bray.csv"))
write.csv(as.data.frame(perm_jac),  file.path(out.dir, "table_permanova_jaccard.csv"))

# 5-3. PERMDISP
disp_bray <- betadisper(bray, site)
disp_jac  <- betadisper(jac,  site)
permutest(disp_bray, permutations = 999)
permutest(disp_jac,  permutations = 999)

disp_df <- bind_rows(
  data.frame(Distance = disp_bray$distances, Site = site, Index = "Bray-Curtis"),
  data.frame(Distance = disp_jac$distances,  Site = site, Index = "Jaccard"))
