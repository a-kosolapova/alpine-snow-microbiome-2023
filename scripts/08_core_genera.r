# ==============================================================================
# Script 08: CORE GENUS-ENVIRONMENT CORRELATIONS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(Hmisc)
library(patchwork)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Helper function: Flatten correlation matrix ----
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = cormat[ut],
    p = pmat[ut]
  )
}

# Aggregate to genus level and transform to relative abundance ----
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_rel <- transform(ps_genus, transform = "compositional")

# Extract metadata ----
metadata <- data.frame(sample_data(ps_rel))

# Select environmental variables ----
chem_data <- metadata %>%
  select(elevation, pH, sodium, potassium, magnesium, calcium, 
         fluoride, chloride, bromide, sulfate, phosphate, acids, TIN)

# Chemical notation labels ----
label_conversion <- c(
  "sodium" = "Na^'+'",
  "potassium" = "K^'+'",
  "magnesium" = "Mg^'2+'",
  "calcium" = "Ca^'2+'",
  "fluoride" = "F^'-'",
  "chloride" = "Cl^'-'",
  "bromide" = "Br^'-'",
  "sulfate" = "SO[4]^'2-'",
  "phosphate" = "PO[4]^'3-'",
  "acids" = "'Organic~Acids'",
  "TIN" = "'TIN'",
  "elevation" = "Elevation",
  "pH" = "pH"
)


# Surface layer -----------------------------------------------------------

# Subset Surface samples
ps_surface <- subset_samples(ps_rel, layer == "Surface")

# Identify core genera (detection = 0.005, prevalence = 50%)
ps_core_surface <- core(ps_surface, detection = 0.005, prevalence = 0.5)

# Extract core genus abundances
tax_surface <- as.data.frame(tax_table(ps_core_surface))
otu_surface <- as.data.frame(otu_table(ps_core_surface))

# Prepare data: rows = samples, columns = genera
core_abund_surface <- t(otu_surface)
colnames(core_abund_surface) <- tax_surface$Genus

# Get metadata for Surface samples
chem_surface <- chem_data[rownames(core_abund_surface), ]

# Combine chemical data and genus abundances
data_surface <- cbind(chem_surface, core_abund_surface)

# Standardize all variables
data_surface_std <- decostand(data_surface, "standardize", 2)

# Calculate Spearman correlations
cor_surface <- cor(data_surface_std, use = "complete.obs")
res_surface <- rcorr(as.matrix(data_surface_std), type = "spearman")
res_surface$Padj <- matrix(p.adjust(res_surface$P, method = "BH"), ncol = ncol(cor_surface))

# Flatten correlation matrix
flat_corr_surface <- flattenCorrMatrix(res_surface$r, res_surface$Padj)
colnames(flat_corr_surface) <- c("Genus", "Env", "r", "Padj")

# Filter for genus-environment pairs only
genus_names_surface <- colnames(core_abund_surface)
env_vars <- colnames(chem_surface)

flat_corr_surface_filt <- flat_corr_surface %>%
  filter(
    (Genus %in% genus_names_surface & Env %in% env_vars) | 
      (Genus %in% env_vars & Env %in% genus_names_surface)
  ) %>%
  mutate(
    # Ensure genus names are in first column
    Genus_new = ifelse(Genus %in% env_vars, Env, Genus),
    Env_new = ifelse(Genus %in% env_vars, Genus, Env)
  ) %>%
  select(-Genus, -Env) %>%
  rename(Genus = Genus_new, Env = Env_new) %>%
  mutate(Significance = case_when(
    Padj < 0.001 ~ "***",
    Padj < 0.01  ~ "**",
    Padj < 0.05  ~ "*",
    TRUE         ~ ""
  ))

# Plot Surface correlations
p_surface <- ggplot(flat_corr_surface_filt, 
                    aes(x = Genus,
                        y = factor(Env, levels = c("elevation", "pH", "phosphate", "sulfate", 
                                                   "bromide", "chloride", "fluoride", "calcium", 
                                                   "magnesium", "potassium", "sodium", "TIN", "acids")), 
                        fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  geom_text(aes(label = Significance), color = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(fill = "Correlation", title = "Surface Layer") +
  scale_y_discrete(labels = function(x) parse(text = label_conversion[x]))


# Bulk layer --------------------------------------------------------------

# Subset Bulk samples
ps_bulk <- subset_samples(ps_rel, layer == "Bulk")

# Identify core genera
ps_core_bulk <- core(ps_bulk, detection = 0.005, prevalence = 0.5)

# Extract core genus abundances
tax_bulk <- as.data.frame(tax_table(ps_core_bulk))
otu_bulk <- as.data.frame(otu_table(ps_core_bulk))

# Prepare data
core_abund_bulk <- t(otu_bulk)
colnames(core_abund_bulk) <- tax_bulk$Genus

# Get metadata for Bulk samples
chem_bulk <- chem_data[rownames(core_abund_bulk), ]

# Combine and standardize
data_bulk <- cbind(chem_bulk, core_abund_bulk)
data_bulk_std <- decostand(data_bulk, "standardize", 2)

# Calculate Spearman correlations
cor_bulk <- cor(data_bulk_std, use = "complete.obs")
res_bulk <- rcorr(as.matrix(data_bulk_std), type = "spearman")
res_bulk$Padj <- matrix(p.adjust(res_bulk$P, method = "BH"), ncol = ncol(cor_bulk))

# Flatten correlation matrix
flat_corr_bulk <- flattenCorrMatrix(res_bulk$r, res_bulk$Padj)
colnames(flat_corr_bulk) <- c("Genus", "Env", "r", "Padj")

# Filter for genus-environment pairs
genus_names_bulk <- colnames(core_abund_bulk)

flat_corr_bulk_filt <- flat_corr_bulk %>%
  filter(
    (Genus %in% genus_names_bulk & Env %in% env_vars) | 
      (Genus %in% env_vars & Env %in% genus_names_bulk)
  ) %>%
  mutate(
    Genus_new = ifelse(Genus %in% env_vars, Env, Genus),
    Env_new = ifelse(Genus %in% env_vars, Genus, Env)
  ) %>%
  select(-Genus, -Env) %>%
  rename(Genus = Genus_new, Env = Env_new) %>%
  mutate(Significance = case_when(
    Padj < 0.001 ~ "***",
    Padj < 0.01  ~ "**",
    Padj < 0.05  ~ "*",
    TRUE         ~ ""
  ))

# Plot Bulk correlations
p_bulk <- ggplot(flat_corr_bulk_filt, 
                 aes(x = Genus,
                     y = factor(Env, levels = c("elevation", "pH", "phosphate", "sulfate", 
                                                "bromide", "chloride", "fluoride", "calcium", 
                                                "magnesium", "potassium", "sodium", "TIN", "acids")), 
                     fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  geom_text(aes(label = Significance), color = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(fill = "Correlation (r)", title = "Bulk Layer") +
  scale_y_discrete(labels = function(x) parse(text = label_conversion[x]))

# Combine plots
combined_plot <- p_surface + p_bulk + plot_layout(widths = c(1.9, 1))

# Save combined plot
ggsave("results/figures/core_genus_correlation_by_layer.png",
       plot = combined_plot,
       height = 6, width = 12, bg = "white", dpi = 600)

# Export correlation tables
write.csv(flat_corr_surface_filt, "results/tables/core_genus_correlations_surface.csv", row.names = FALSE)
write.csv(flat_corr_bulk_filt, "results/tables/core_genus_correlations_bulk.csv", row.names = FALSE)