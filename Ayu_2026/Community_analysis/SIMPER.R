#-------------------------------------------------------------------------------
# 8.SIMPER（どの属が地点間の違いに効いているか）
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(ggplot2)
library(patchwork)

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
# 注意: SIMPERは「違いに効いている属の目安」で、有意差の検定ではない
set.seed(123)
sim <- simper(spe, site, permutations = 99)
sim_sum <- summary(sim, ordered = TRUE)
sim_sum

# 8-3. 地点の組み合わせごとに上位5属を保存
top5 <- do.call(rbind, lapply(names(sim_sum), function(nm) {
  d <- head(sim_sum[[nm]], 5)
  data.frame(Pair = nm, Genus = rownames(d), contribution = d$average, cumsum = d$cumsum)
}))
top5
write.csv(top5, file.path(out.dir, "table_simper_top5.csv"), row.names = FALSE)

# 8-4. Barplot（全ペアを1枚に）
# 1つの地点ペアの上位5属を、寄与率の大きい順に上から並べる関数
plot_simper_pair <- function(d, pair_name) {
  d$Genus <- factor(d$Genus, levels = rev(d$Genus[order(d$contribution, decreasing = TRUE)]))
  ggplot(d, aes(contribution * 100, Genus)) +
    geom_col(fill = "#4C78A8", color = "black", linewidth = 0.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Contribution to dissimilarity (%)", y = NULL, title = pair_name) +
    theme_classic(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

pairs_use <- unique(as.character(top5$Pair))   # 全ペア(4地点なら6つ)
plots <- lapply(pairs_use, function(p) plot_simper_pair(top5[top5$Pair == p, ], p))

n_pair <- length(plots)
fig_simper <- wrap_plots(plots, ncol = 2) +
  plot_annotation(tag_levels = "A")   # 各パネルに A, B, C... を付ける
print(fig_simper)

ggsave(file.path(out.dir, "fig8_simper_top5_all_pairs.png"), fig_simper,
       width = 11, height = ceiling(n_pair / 2) * 4, dpi = 300, bg = "white")