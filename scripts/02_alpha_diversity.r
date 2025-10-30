# ==============================================================================
# Script 02: ALPHA-DIVERSITY ANALYSIS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(rstatix)
library(biomeUtils)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Calculate alpha diversity indices ----
alpha_div <- estimate_richness(ps)

# Merge alpha diversity with metadata ----
metadata_alpha <- cbind(
  as.data.frame(sample_data(ps)),
  alpha_div
)

# Calculate Faith's Phylogenetic Diversity ----
ps_pd <- calculatePD(ps)
pd_data <- as.data.frame(sample_data(ps_pd))
metadata_alpha$PD <- pd_data$PD

# Convert layer to factor ----
metadata_alpha$layer <- as.factor(metadata_alpha$layer)

# Color palette for layers ----
layer_colors <- c("Surface" = "blue", "Bulk" = "orange")

# Shannon Diversity ----
# Statistical test
shannon_test <- pairwise_wilcox_test(metadata_alpha, Shannon ~ layer)
shannon_test <- add_xy_position(shannon_test, x = "layer")

# Plot
p_shannon <- ggplot(metadata_alpha, aes(x = layer, y = Shannon)) +  
  geom_boxplot(aes(fill = layer), alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  labs(fill = "Snow layer") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("Snow layer") +
  ylab("Shannon Diversity") +
  scale_fill_manual(values = layer_colors) +
  stat_pvalue_manual(shannon_test, 
                     label = "p.adj.signif",
                     y.position = max(metadata_alpha$Shannon) * 1.1)

ggsave("results/figures/alpha_diversity_shannon.png", 
       plot = p_shannon, 
       height = 6, width = 3, bg = "white", dpi = 600)

# Inverse Simpson Diversity ----
# Statistical test
invsimp_test <- pairwise_wilcox_test(metadata_alpha, InvSimpson ~ layer)
invsimp_test <- add_xy_position(invsimp_test, x = "layer")

# Plot
p_invsimp <- ggplot(metadata_alpha, aes(x = layer, y = InvSimpson)) +  
  geom_boxplot(aes(fill = layer), alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  labs(fill = "Snow layer") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("Snow layer") +
  ylab("Inverse Simpson Diversity") +
  scale_fill_manual(values = layer_colors) +
  stat_pvalue_manual(invsimp_test, 
                     label = "p.adj.signif",
                     y.position = max(metadata_alpha$InvSimpson) * 1.1)

ggsave("results/figures/alpha_diversity_invsimpson.png", 
       plot = p_invsimp, 
       height = 6, width = 3, bg = "white", dpi = 600)

# Faith's Phylogenetic Diversity ----
# Statistical test
pd_test <- pairwise_wilcox_test(metadata_alpha, PD ~ layer)
pd_test <- add_xy_position(pd_test, x = "layer")

# Plot
p_pd <- ggplot(metadata_alpha, aes(x = layer, y = PD)) +  
  geom_boxplot(aes(fill = layer), alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  labs(fill = "Snow layer") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("Snow layer") +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = layer_colors) 

ggsave("results/figures/alpha_diversity_faiths_pd.png", 
       plot = p_pd, 
       height = 6, width = 3, bg = "white", dpi = 600)



