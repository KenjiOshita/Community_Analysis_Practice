#-------------------------------------------------------------------------------
# 4.NMDS (Bray-Curtis and Jaccard)
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "NMDS"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 4-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)
site_cols  <- setNames(brewer.pal(8, "Set2")[(seq_along(site_order) - 1) %% 8 + 1], site_order)

# 4-2. NMDS
# Bray-Curtis: 量も見る / Jaccard: いる・いないだけ見る
set.seed(123)   # NMDSは乱数を使うので固定
nmds_bray <- metaMDS(spe, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
nmds_jac  <- metaMDS(spe, distance = "jaccard", binary = TRUE, k = 2, trymax = 100,
                     autotransform = FALSE)
nmds_bray$stress   # 0.2未満ならOK、0.1未満なら良好
nmds_jac$stress

# 4-3. Shepard plot
png(file.path(out.dir, "fig4a_nmds_shepard.png"), width = 2800, height = 1400, res = 300)
par(mfrow = c(1, 2))
stressplot(nmds_bray, main = "Bray-Curtis")
stressplot(nmds_jac,  main = "Jaccard")
dev.off()

# 4-4. Plot
make_df <- function(nmds, label) {
  sc <- scores(nmds, display = "sites")
  data.frame(NMDS1 = sc[, 1], NMDS2 = sc[, 2], SampleID = rownames(sc),
             Site = factor(sub("_[0-9]+$", "", rownames(sc)), levels = site_order),
             Index = paste0(label, " (stress = ", round(nmds$stress, 3), ")"))
}
nmds_df <- bind_rows(make_df(nmds_bray, "Bray-Curtis"),
                     make_df(nmds_jac,  "Jaccard"))

fig_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, fill = Site)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  stat_ellipse(aes(color = Site), geom = "polygon", alpha = 0.12, linewidth = 0.3) +
  geom_point(shape = 21, size = 3.5, color = "black") +
  facet_wrap(~ Index) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  coord_fixed() +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = "black"))
print(fig_nmds)
ggsave(file.path(out.dir, "fig4b_nmds.png"), fig_nmds,
       width = 11, height = 5.5, dpi = 300, bg = "white")
