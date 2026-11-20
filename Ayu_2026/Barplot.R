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
