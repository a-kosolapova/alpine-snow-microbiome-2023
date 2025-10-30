# ==============================================================================
# Script 04: HIERARCHICAL CLUSTERING AND TAXONOMIC COMPOSITION

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(vegan)
library(microbiome)
library(ggdendro)
library(cowplot)
library(colorspace)
library(cluster)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Calculate hierarchical clustering ----
# Hellinger transform
ps_hell <- transform(ps, transform = "hellinger")

# Calculate weighted UniFrac distance
dist_wunifrac <- UniFrac(ps_hell, weighted = TRUE, normalized = TRUE, 
                         parallel = FALSE, fast = TRUE)

# Perform hierarchical clustering using agnes with flexible beta method
cluster_flexible <- agnes(x = dist_wunifrac, method = "flexible", par.method = 0.625)
cluster_flexible_hclust <- as.hclust(cluster_flexible)

# Create dendrogram object
dend <- as.dendrogram(cluster_flexible_hclust)

# Taxonomic level for visualization ----
tax_level <- "Phylum" 

# Transform to relative abundance ----
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Aggregate to selected taxonomic level ----
ps_tax <- tax_glom(ps_rel, taxrank = tax_level)

# Melt data to long format ----
df_long <- psmelt(ps_tax)

# Define core taxa ----
core_taxa <- core(ps_rel, detection = 0.002, prevalence = 0.25)
core_taxa_names <- as.character(tax_table(core_taxa)[, tax_level])
core_taxa_names <- unique(core_taxa_names[!is.na(core_taxa_names)])

# Collapse non-core taxa into "Other" ----
df_long <- df_long %>%
  mutate(!!tax_level := ifelse(!!sym(tax_level) %in% core_taxa_names, 
                               as.character(!!sym(tax_level)), "Other"))

# Re-aggregate after collapsing ----
df_wide <- df_long %>%
  group_by(Sample, !!sym(tax_level)) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = !!sym(tax_level), values_from = Abundance, values_fill = 0)

# Order samples by dendrogram ----
sample_order <- labels(dend)
df_wide <- df_wide %>% 
  column_to_rownames("Sample") %>%
  .[sample_order, ]

# Prepare plotting data ----
dend_data <- dendro_data(dend, type = "rectangle")
segment_data <- dend_data[["segments"]]
sample_pos_table <- with(dend_data$labels,
                         data.frame(x_center = x, Sample = as.character(label), width = 0.9))

# Prepare stacked bar data ----
df_long_stacked <- df_wide %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = tax_level, values_to = "Abundance") %>%
  group_by(Sample) %>%
  mutate(
    Frequency = Abundance * 100,
    ymax = cumsum(Frequency / sum(Frequency)),
    ymin = lag(ymax, default = 0),
    y_center = (ymax + ymin) / 2
  ) %>%
  ungroup() %>%
  left_join(sample_pos_table, by = "Sample") %>%
  mutate(
    xmin = x_center - width / 2,
    xmax = x_center + width / 2
  )

# Set axis limits ----
buffer <- 0.5
axis_limits <- with(sample_pos_table,
                    c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))) +
  buffer * c(-1, 1)

# Prepare sample metadata for annotation ----
sample_metadata <- data.frame(
  Sample = sample_names(ps),
  Layer = as(sample_data(ps)$layer, "character"),
  Valley = as(sample_data(ps)$valley, "character")
)

sample_metadata <- sample_metadata %>%
  mutate(
    Valley_simple = case_when(
      Valley == "Val de Valsorrey" ~ "Valley1",
      Valley == "Val de Bagnes" ~ "Valley2",
      Valley == "Val Ferret" ~ "Valley3",
      TRUE ~ "Valley4"
    ),
    Layer_Valley = paste(Layer, Valley_simple, sep = "_")
  )

# Define color palette ----
base_colors <- c("Surface" = "#1f77b4", "Bulk" = "#ff7f0e")

layer_valley_colors <- c(
  "Surface_Valley1" = lighten(unname(base_colors["Surface"]), amount = 0.3),
  "Surface_Valley2" = lighten(unname(base_colors["Surface"]), amount = 0.15),
  "Surface_Valley3" = darken(unname(base_colors["Surface"]), amount = 0.15),
  "Surface_Valley4" = darken(unname(base_colors["Surface"]), amount = 0.3),
  "Bulk_Valley1" = lighten(unname(base_colors["Bulk"]), amount = 0.3),
  "Bulk_Valley2" = lighten(unname(base_colors["Bulk"]), amount = 0.15),
  "Bulk_Valley3" = darken(unname(base_colors["Bulk"]), amount = 0.15),
  "Bulk_Valley4" = darken(unname(base_colors["Bulk"]), amount = 0.3)
)

# Create layer annotation ----
layer_bar_df <- left_join(sample_pos_table, sample_metadata, by = "Sample")

layer_annotation <- ggplot(layer_bar_df, aes(x = x_center, fill = Layer_Valley)) +
  geom_tile(aes(y = 0, height = 0.1, width = width), color = "black") +
  scale_fill_manual(values = layer_valley_colors) +
  scale_x_continuous(limits = axis_limits, expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
    legend.position = "none"
  )

# Create dendrogram plot ----
plt_dendro <- ggplot(segment_data) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_y_continuous(expand = c(0, 0.05)) +
  scale_x_continuous(breaks = sample_pos_table$x_center,
                     labels = rep("", nrow(sample_pos_table)),
                     limits = axis_limits,
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Distance") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Create relative abundance plot ----
plt_abundance <- ggplot(df_long_stacked,
                        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                            fill = .data[[tax_level]])) +
  geom_rect(color = "white", linewidth = 0.2) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_x_continuous(breaks = sample_pos_table$x_center,
                     labels = sample_pos_table$Sample,
                     limits = axis_limits,
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Relative Abundance", fill = tax_level) +
  theme_bw() +
  theme(
    plot.margin = unit(c(-0.9, 0.2, 1, 0.2), "cm"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.box = "horizontal"
  ) +
  scale_fill_brewer(palette = "Set3")

# Combine plots ----
final_plot <- plot_grid(
  plt_dendro,
  layer_annotation,
  plt_abundance,
  align = "v",
  ncol = 1,
  rel_heights = c(0.4, 0.05, 1)
)

# Save plot ----
ggsave("results/figures/fig3.png",
       plot = final_plot,
       height = 10, width = 14, bg = "white", dpi = 600)

