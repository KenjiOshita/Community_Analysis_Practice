#-------------------------------------------------------------------------------
# 6.Alpha diversity (Richness and Shannon)
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "Alpha"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 6-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)   # 上流→下流の順に並べたいとき: c("Upper","Middle","Lower","Estuary")
site_cols  <- setNames(brewer.pal(8, "Set2")[(seq_along(site_order) - 1) %% 8 + 1], site_order)

# 6-2. Calculate indices
alpha_df <- data.frame(
  SampleID = rownames(spe),
  Site     = factor(site_names, levels = site_order),
  Richness = specnumber(spe),                  # 属数
  Shannon  = diversity(spe, index = "shannon") # シャノン指数 H'
)
write.csv(alpha_df, file.path(out.dir, "table_alpha_diversity.csv"), row.names = FALSE)

# 6-3. Plot: 箱ひげ + 個体の点 + 地点平均の折れ線（白いひし形）
alpha_long <- alpha_df %>%
  pivot_longer(c(Richness, Shannon), names_to = "Index", values_to = "Value")
alpha_long$Index <- factor(alpha_long$Index, levels = c("Richness", "Shannon"),
                           labels = c("Genus richness", "Shannon index (H')"))

alpha_mean <- alpha_long %>%
  group_by(Index, Site) %>%
  summarise(Mean = mean(Value), .groups = "drop")

fig_alpha <- ggplot(alpha_long, aes(Site, Value)) +
  geom_boxplot(aes(fill = Site), alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, color = "grey30") +
  geom_line(data = alpha_mean, aes(Site, Mean, group = 1), linewidth = 0.9) +
  geom_point(data = alpha_mean, aes(Site, Mean), shape = 23, size = 3.5, fill = "white") +
  facet_wrap(~ Index, scales = "free_y") +
  scale_fill_manual(values = site_cols) +
  labs(x = "Site", y = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = "black"))
print(fig_alpha)
ggsave(file.path(out.dir, "fig6_alpha_diversity.png"), fig_alpha,
       width = 9, height = 4.5, dpi = 300, bg = "white")

# 6-4. Statistics: Kruskal-Wallis -> pairwise Wilcoxon (Holm)
kruskal.test(Richness ~ Site, data = alpha_df)
kruskal.test(Shannon  ~ Site, data = alpha_df)
pairwise.wilcox.test(alpha_df$Richness, alpha_df$Site, p.adjust.method = "holm", exact = FALSE)
pairwise.wilcox.test(alpha_df$Shannon,  alpha_df$Site, p.adjust.method = "holm", exact = FALSE)

