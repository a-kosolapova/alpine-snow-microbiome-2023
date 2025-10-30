# ==============================================================================
# SCRIPT 05: DISTANCE-BASED REDUNDANCY ANALYSIS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(vegan)
library(microbiome)
library(patchwork)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Extract metadata ----
metadata <- data.frame(sample_data(ps))

# Hellinger transform ----
ps_hell <- transform(ps, transform = "hellinger")

# Extract Hellinger-transformed OTU table ----
otu_hell <- t(as.matrix(otu_table(ps_hell)))

# Calculate distances ----
# Bray-Curtis on Hellinger-transformed data
bc_dist_hel <- vegdist(otu_hell, method = "bray")

# Weighted UniFrac
wunifrac_dist <- UniFrac(ps_hell, weighted = TRUE, normalized = TRUE, 
                         parallel = FALSE, fast = TRUE)

# Log-transform chemical variables ----
log_transform <- function(x) {
  log(x + (min(x[x > 0]) / 2))
}

metadata$logPotassium <- log_transform(metadata$potassium)
metadata$logSodium <- log_transform(metadata$sodium)
metadata$logMagnesium <- log_transform(metadata$magnesium)
metadata$logPhosphate <- log_transform(metadata$phosphate)
metadata$logSulfate <- log_transform(metadata$sulfate)
metadata$logCalcium <- log_transform(metadata$calcium)
metadata$logAcids <- log_transform(metadata$acids)
metadata$logTIN <- log_transform(metadata$TIN)
metadata$logChloride <- log_transform(metadata$chloride)
metadata$logElevation <- log_transform(metadata$elevation)

# Select environmental variables ----
selected_chemicals <- metadata %>%
  select(logPotassium, logSodium, logMagnesium, logPhosphate,logElevation,
         logSulfate, logCalcium, logAcids, logTIN, logChloride, pH)


# Bray-Curtis dbRDA -------------------------------------------------------

# Define null and full models for Bray-Curtis
null_model_bc <- capscale(bc_dist_hel ~ 1, data = selected_chemicals)
full_model_bc <- capscale(bc_dist_hel ~ ., data = selected_chemicals)

# Model selection with forward/backward stepwise
set.seed(123)
selected_model_bc <- ordistep(null_model_bc, 
                              scope = formula(full_model_bc), 
                              direction = "both", 
                              permutations = 999)

# Test model significance
anova(selected_model_bc)
anova.cca(selected_model_bc, by = "margin")

# Extract variance components
total_inertia_bc <- selected_model_bc$tot.chi
constrained_inertia_bc <- sum(selected_model_bc$CCA$eig)
eigenvals_bc <- selected_model_bc$CCA$eig

# Percent of fitted and total variation
var_fitted_bc <- round(100 * eigenvals_bc / constrained_inertia_bc, 1)
var_total_bc <- round(100 * eigenvals_bc / total_inertia_bc, 1)

# Extract site scores
site_scores_bc <- scores(selected_model_bc, display = "sites")
dbrda_points_bc <- as.data.frame(site_scores_bc)
dbrda_points_bc$valley <- metadata$valley
dbrda_points_bc$layer <- metadata$layer

# Extract environmental vectors
env_vectors_bc <- scores(selected_model_bc, display = "bp")
env_df_bc <- as.data.frame(env_vectors_bc)


# wUniFrac dbRDA ----------------------------------------------------------

# Define null and full models for weighted UniFrac
null_model_wuf <- capscale(wunifrac_dist ~ 1, data = selected_chemicals)
full_model_wuf <- capscale(wunifrac_dist ~ ., data = selected_chemicals)

# Model selection with forward/backward stepwise
set.seed(123)
selected_model_wuf <- ordistep(null_model_wuf, 
                               scope = formula(full_model_wuf), 
                               direction = "both", 
                               permutations = 999)

# Test model significance
anova(selected_model_wuf)
anova.cca(selected_model_wuf, by = "margin")

# Extract variance components
total_inertia_wuf <- selected_model_wuf$tot.chi
constrained_inertia_wuf <- sum(selected_model_wuf$CCA$eig)
eigenvals_wuf <- selected_model_wuf$CCA$eig

# Percent of fitted and total variation
var_fitted_wuf <- round(100 * eigenvals_wuf / constrained_inertia_wuf, 1)
var_total_wuf <- round(100 * eigenvals_wuf / total_inertia_wuf, 1)

# Extract site scores
site_scores_wuf <- scores(selected_model_wuf, display = "sites")
dbrda_points_wuf <- as.data.frame(site_scores_wuf)
dbrda_points_wuf$valley <- metadata$valley
dbrda_points_wuf$layer <- metadata$layer

# Extract environmental vectors
env_vectors_wuf <- scores(selected_model_wuf, display = "bp")
env_df_wuf <- as.data.frame(env_vectors_wuf)


# Prepare data for plotting -----------------------------------------------

# Chemical notation labels
label_conversion <- c(
  "logSodium" = "Na^'+'",
  "logPotassium" = "K^'+'",
  "logMagnesium" = "Mg^'2+'",
  "logCalcium" = "Ca^'2+'",
  "logChloride" = "Cl^'-'",
  "logSulfate" = "SO[4]^'2-'",
  "logPhosphate" = "PO[4]^'3-'",
  "logAcids" = "'Organic\\nAcids'",
  "logTIN" = "'TIN'",
  "pH" = "pH",
  "logElevation" = "'Elevation'"
)

# Apply labels to Bray-Curtis
pretty_names_bc <- label_conversion[rownames(env_df_bc)]
pretty_names_bc[is.na(pretty_names_bc)] <- rownames(env_df_bc)[is.na(pretty_names_bc)]
env_df_bc$var <- pretty_names_bc

# Apply labels to weighted UniFrac
pretty_names_wuf <- label_conversion[rownames(env_df_wuf)]
pretty_names_wuf[is.na(pretty_names_wuf)] <- rownames(env_df_wuf)[is.na(pretty_names_wuf)]
env_df_wuf$var <- pretty_names_wuf

# Arrow scaling
arrow_scale <- 3


# Create plots ------------------------------------------------------------

# Bray-Curtis dbRDA plot
BC_dbRDA <- ggplot(dbrda_points_bc, aes(x = CAP1, y = CAP2, color = layer)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_segment(
    data = env_df_bc,
    aes(x = 0, y = 0, xend = CAP1 * arrow_scale, yend = CAP2 * arrow_scale),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = env_df_bc,
    aes(x = CAP1 * arrow_scale * 1.1,
        y = CAP2 * arrow_scale * 1.1, 
        label = var),
    parse = TRUE,  
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = paste0("dbRDA1 (", var_fitted_bc[1], "% of fitted, ", var_total_bc[1], "% of total)"),
    y = paste0("dbRDA2 (", var_fitted_bc[2], "% of fitted, ", var_total_bc[2], "% of total)"),
    color = "Layer",
    title = "Bray-Curtis"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Surface" = "blue", "Bulk" = "orange")) +
  stat_ellipse(aes(group = layer), linetype = "dashed", level = 0.95) + 
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "none"
  )

# Weighted UniFrac dbRDA plot
wUniFrac_dbRDA <- ggplot(dbrda_points_wuf, aes(x = CAP1, y = CAP2, color = layer)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_segment(
    data = env_df_wuf,
    aes(x = 0, y = 0, xend = CAP1 * arrow_scale, yend = CAP2 * arrow_scale),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = env_df_wuf,
    aes(x = CAP1 * arrow_scale * 1.1,
        y = CAP2 * arrow_scale * 1.1, 
        label = var),
    parse = TRUE,  
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = paste0("dbRDA1 (", var_fitted_wuf[1], "% of fitted, ", var_total_wuf[1], "% of total)"),
    y = paste0("dbRDA2 (", var_fitted_wuf[2], "% of fitted, ", var_total_wuf[2], "% of total)"),
    color = "Layer",
    title = "Weighted UniFrac"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Surface" = "blue", "Bulk" = "orange")) +
  stat_ellipse(aes(group = layer), linetype = "dashed", level = 0.95) + 
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "none"
  )

# Combine plots
combined_plot <-  wUniFrac_dbRDA + BC_dbRDA 

# Save combined plot
ggsave("results/figures/dbrda_combined.png", 
       plot = combined_plot, 
       height = 6, width = 12, bg = "white", dpi = 600)

