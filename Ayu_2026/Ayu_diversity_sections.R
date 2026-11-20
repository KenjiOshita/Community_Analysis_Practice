#Title Alpha / Beta diversity of algae in ayu gut contents (genus level)
#Edit day 2026/11/20
#
# 【このファイルの使い方】
#  ・各セクション（1〜9）は「独立」しています。どのセクションも
#    rm() → library → ファイル読込 から始まるので、必要なセクションだけ選んで実行できます。
#  ・ただし 2〜9 は、セクション1が作る「Ayu_Rarefied_Table.csv」を読み込みます。
#    → まず必ずセクション1を1回実行すること。
#  ・使うファイルは2つ（CSV）:
#       Ayu_counts.csv : 1列目 SeqName、2列目以降 = サンプルごとのリード数
#       Ayu_taxo.csv   : SeqName と Genus を含む分類表
#  ・サンプル名は「地点名_番号」（例 Upper_01）にしておく。
#    地点は、サンプル名の最後の「_番号」を除いた部分から自動で取り出す。
#      セクション一覧
#       1. レアファクション      6. α多様性
#       2. 組成グラフ（属）      7. ペアワイズPERMANOVA
#       3. クラスタリング        8. SIMPER
#       4. NMDS                  9. envfit（NMDS + 属の矢印）
#       5. PERMANOVA・PERMDISP
#  ・必要なパッケージは、下の「0.Install packages」で、入っていないものだけ自動で入る。
#-------------------------------------------------------------------------------


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
site_cols  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

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


#-------------------------------------------------------------------------------
# 2.Genus composition bar plot
#-------------------------------------------------------------------------------
rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "Barplot"
if (!dir.exists(out.dir)) dir.create(out.dir)

min_pct <- 3   # どこかで この%以上 にならない属は Others にまとめる

# 2-1. Import data and convert to relative abundance (%)
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)

genus_long <- raw %>%
  pivot_longer(cols = -SampleID, names_to = "Genus", values_to = "count") %>%
  group_by(SampleID) %>%
  mutate(Percent = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(Site = sub("_[0-9]+$", "", SampleID))
head(genus_long)

site_order <- unique(genus_long$Site)
genus_long$Site <- factor(genus_long$Site, levels = site_order)

# 2-2. Others にまとめる
keep_genus <- genus_long %>%
  group_by(Genus) %>%
  summarise(max_pct = max(Percent)) %>%
  filter(max_pct >= min_pct) %>%
  pull(Genus)
keep_genus <- sort(keep_genus)

genus_long <- genus_long %>%
  mutate(Genus = ifelse(Genus %in% keep_genus, Genus, "Others")) %>%
  group_by(SampleID, Site, Genus) %>%
  summarise(Percent = sum(Percent), .groups = "drop")

genus_long$Genus <- factor(genus_long$Genus, levels = c(keep_genus, "Others"))

# 2-3. Colors（Others は灰色）
genus_cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(keep_genus))
names(genus_cols) <- keep_genus
genus_cols["Others"] <- "grey70"

# 2-4. Plot A: individual fish
fig_ind <- ggplot(genus_long, aes(SampleID, Percent, fill = Genus)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ Site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_cols) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
  labs(x = "Individual fish", y = "Relative abundance (%)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey92"))
print(fig_ind)
ggsave(file.path(out.dir, "fig2a_barplot_individual.png"), fig_ind,
       width = 12, height = 6, dpi = 300, bg = "white")

# 2-5. Plot B: site mean
genus_site <- genus_long %>%
  group_by(Site, Genus) %>%
  summarise(Percent = mean(Percent), .groups = "drop")

fig_site <- ggplot(genus_site, aes(Site, Percent, fill = Genus)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = genus_cols) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
  labs(x = "Site", y = "Mean relative abundance (%)") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(face = "italic"),
        panel.grid = element_blank())
print(fig_site)
ggsave(file.path(out.dir, "fig2b_barplot_site_mean.png"), fig_site,
       width = 7, height = 6, dpi = 300, bg = "white")


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
site_cols  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

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
site_cols  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

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
site_cols <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

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

fig_disp <- ggplot(disp_df, aes(Site, Distance, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  facet_wrap(~ Index) +
  scale_fill_manual(values = site_cols) +
  labs(x = "Site", y = "Distance to site centroid") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = "black"))
print(fig_disp)
ggsave(file.path(out.dir, "fig5_permdisp.png"), fig_disp,
       width = 9, height = 4.5, dpi = 300, bg = "white")


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
site_cols  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

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


#-------------------------------------------------------------------------------
# 9.envfit（NMDSに、効いている属を矢印で重ねる）
#-------------------------------------------------------------------------------
rm(list = ls())

library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

Sys.setenv(LANG = "en")

out.dir <- "NMDS"
if (!dir.exists(out.dir)) dir.create(out.dir)

p_cut <- 0.01   # 矢印を描く属のp値の基準

# 9-1. Import data
raw <- read.csv("Ayu_Rarefied_Table.csv", check.names = FALSE)
spe <- raw[, -1]
rownames(spe) <- raw$SampleID

site_names <- sub("_[0-9]+$", "", rownames(spe))
site_order <- unique(site_names)
site_cols  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(site_order)), site_order)

# 9-2. NMDS (Bray-Curtis) and envfit
set.seed(123)
nmds <- metaMDS(spe, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
ef   <- envfit(nmds, spe, permutations = 999)

arrow_df <- as.data.frame(scores(ef, display = "vectors"))
arrow_df$Genus <- rownames(arrow_df)
arrow_df$p     <- ef$vectors$pvals
arrow_df       <- filter(arrow_df, p < p_cut)
arrow_df

sc <- scores(nmds, display = "sites")
site_df <- data.frame(NMDS1 = sc[, 1], NMDS2 = sc[, 2],
                      Site = factor(site_names, levels = site_order))

# 9-3. Plot
fig_envfit <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point(data = site_df, aes(NMDS1, NMDS2, fill = Site),
             shape = 21, size = 3.5, color = "black") +
  geom_segment(data = arrow_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "grey30") +
  geom_text_repel(data = arrow_df, aes(NMDS1, NMDS2, label = Genus),
                  fontface = "italic", size = 4, segment.color = NA) +
  scale_fill_manual(values = site_cols) +
  coord_fixed() +
  labs(caption = paste("Bray-Curtis NMDS, stress =", round(nmds$stress, 3))) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold"))
print(fig_envfit)
ggsave(file.path(out.dir, "fig7_nmds_envfit.png"), fig_envfit,
       width = 7, height = 6, dpi = 300, bg = "white")
