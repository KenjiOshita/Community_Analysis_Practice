#-------------------------------------------------------------------------------
# 0.Install packages（最初に1回だけ実行。入っていないパッケージだけをインストール）
#-------------------------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "vegan",
          "RColorBrewer", "dendextend", "ggrepel")

new_pkgs <- pkgs[!(pkgs %in% rownames(installed.packages()))]   # 未インストールのもの
if (length(new_pkgs) > 0) install.packages(new_pkgs)


#-------------------------------------------------------------------------------
# 1.Rarefaction
#-------------------------------------------------------------------------------
rm(list = ls())

library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "Rarefaction"
if (!dir.exists(out.dir)) dir.create(out.dir)

# 1-1. Import data
count_data <- read.csv("Ayu_counts.csv", check.names = FALSE)
taxa_data  <- read.csv("Ayu_taxo.csv")
head(count_data)
head(taxa_data)

# 1-2. ASV -> Genus table
# 属が空欄のASVは除き、SeqNameで分類とリード数をくっつけて、同じ属を足し算する
taxa_use <- taxa_data %>%
  filter(!is.na(Genus), Genus != "") %>%
  select(SeqName, Genus)

genus_table <- inner_join(taxa_use, count_data, by = "SeqName") %>%
  select(-SeqName) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

# 行 = サンプル、列 = 属 にひっくり返す（vegan用の形）
genus_mat <- genus_table %>%
  column_to_rownames("Genus") %>%
  as.matrix() %>%
  t()
head(genus_mat)

# 1-3. リード数の確認
depth <- rowSums(genus_mat)
depth
min_depth <- min(depth)   # 最も少ないサンプルに合わせる
min_depth

# 1-4. Rarefaction curve
rare_curve <- rarecurve(genus_mat, step = 200, tidy = TRUE)
names(rare_curve) <- c("SampleID", "Reads", "Richness")
rare_curve$Site <- sub("_[0-9]+$", "", rare_curve$SampleID)

site_order <- unique(rare_curve$Site)
site_cols  <- setNames(brewer.pal(8, "Set2")[(seq_along(site_order) - 1) %% 8 + 1], site_order)

fig_rare <- ggplot(rare_curve, aes(Reads, Richness, group = SampleID, color = Site)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = min_depth, linetype = "dashed") +
  scale_color_manual(values = site_cols) +
  labs(x = "Number of reads", y = "Number of genera",
       caption = paste("Dashed line: rarefaction depth =", min_depth)) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"))
print(fig_rare)
ggsave(file.path(out.dir, "fig1_rarefaction_curve.png"), fig_rare,
       width = 7, height = 5, dpi = 300, bg = "white")

# 1-5. Rarefy to the minimum depth
# ランダムに抜き取るので、乱数の種を固定しておく（毎回同じ結果になる）
set.seed(123)
rare_mat <- rrarefy(genus_mat, sample = min_depth)
rare_mat <- rare_mat[, colSums(rare_mat) > 0]   # 全部0になった属は除く
rowSums(rare_mat)                               # 全て同じ数ならOK

# 1-6. Save (以降のセクションはこのファイルを読み込む)
output <- data.frame(SampleID = rownames(rare_mat), rare_mat, check.names = FALSE)
write.csv(output, "Ayu_Rarefied_Table.csv", row.names = FALSE)

# 1-7. レアファクションカーブ：前と後の比較（任意。論文・スライド用）
# 前（Before）: 元のリード数のまま
# 後（After） : 最小リード数にそろえた後。全サンプルが同じ深さで止まる
rare_curve_after <- rarecurve(rare_mat, step = 200, tidy = TRUE)
names(rare_curve_after) <- c("SampleID", "Reads", "Richness")
rare_curve_after$Site <- sub("_[0-9]+$", "", rare_curve_after$SampleID)

# 2枚の図で軸をそろえる（そろえないと見比べられない）
x_max <- max(rare_curve$Reads)
y_max <- max(rare_curve$Richness)

fig_rare_before <- ggplot(rare_curve, aes(Reads, Richness, group = SampleID, color = Site)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = site_cols) +
  coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
  labs(title = "Before rarefaction", x = "Number of reads", y = "Number of genera") +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"))

fig_rare_after <- ggplot(rare_curve_after, aes(Reads, Richness, group = SampleID, color = Site)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = site_cols) +
  coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
  labs(title = "After rarefaction", x = "Number of reads", y = "Number of genera") +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"))

print(fig_rare_before)
print(fig_rare_after)
ggsave(file.path(out.dir, "fig1b_rarefaction_curve_before.png"), fig_rare_before,
       width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out.dir, "fig1c_rarefaction_curve_after.png"), fig_rare_after,
       width = 7, height = 5, dpi = 300, bg = "white")

