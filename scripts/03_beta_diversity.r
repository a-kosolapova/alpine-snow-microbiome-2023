# ============================================================================
# Script 04: BETA-DIVERSITY ANALYSIS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(rstatix)
library(microbiome)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Hellinger transform (use for both Bray-Curtis and UniFrac) ----
ps_hell <- transform(ps, transform = "hellinger")

# Extract OTU table and metadata ----
otu_table_hel <- t(as.matrix(otu_table(ps_hell)))
metadata <- data.frame(sample_data(ps_hell))

# Bray-Curtis dissimilarity -----------------------------------------------

# Calculate Bray-Curtis dissimilarity ----
dist_bc <- vegdist(otu_table_hel, method = "bray")

# PERMANOVA test ----
set.seed(4104)
permanova_result <- adonis2(dist_bc ~ layer, data = metadata, permutations = 999)

# Test for homogeneity of dispersions ----
dispersion <- betadisper(dist_bc, metadata$layer)
dispersion_test <- anova(dispersion)

# Convert distance matrix to dataframe ----
dist_bc_matrix <- as.matrix(dist_bc)
dist_bc_df <- data.frame(
  Sample1 = rep(rownames(dist_bc_matrix), each = ncol(dist_bc_matrix)),
  Sample2 = rep(colnames(dist_bc_matrix), times = nrow(dist_bc_matrix)),
  BC = as.vector(dist_bc_matrix)
)

# Remove self-comparisons
dist_bc_df <- dist_bc_df %>% filter(Sample1 != Sample2)

# Function to extract site information
extract_site <- function(sample) {
  sub("_.*", "", sample)
}

# Function to extract layer information
extract_layer <- function(sample) {
  ifelse(grepl("_T$", sample), "Surface",
         ifelse(grepl("_M$", sample), "Bulk", NA))
}

# Add comparison type columns ----
dist_bc_df <- dist_bc_df %>%
  mutate(
    # Determine if comparison is between same layers
    Layer1 = extract_layer(Sample1),
    Layer2 = extract_layer(Sample2),
    ComparisonType = case_when(
      Layer1 == "Surface" & Layer2 == "Surface" ~ "Surface-Surface",
      Layer1 == "Bulk" & Layer2 == "Bulk" ~ "Bulk-Bulk",
      TRUE ~ "Mixed"
    ),
    # Determine if comparison is within or between sites
    Site1 = extract_site(Sample1),
    Site2 = extract_site(Sample2),
    ComparisonType2 = ifelse(Site1 == Site2, "Intra-site", "Inter-site")
  )

# Filter for inter-site comparisons only ----
dist_inter <- dist_bc_df %>%
  filter(ComparisonType2 == "Inter-site") %>%
  filter(ComparisonType != "Mixed")

# Set factor levels
dist_inter$ComparisonType <- factor(dist_inter$ComparisonType, 
                                    levels = c("Surface-Surface", "Bulk-Bulk"))

# Statistical test: Wilcoxon test between Surface and Bulk inter-site distances ----
wilcox_result <- dist_inter %>%
  pairwise_wilcox_test(BC ~ ComparisonType) %>%
  add_xy_position(x = "ComparisonType")


# Weighted UniFrac --------------------------------------------------------

# Calculate weighted UniFrac distance (using same Hellinger-transformed object) ----
dist_wunifrac <- UniFrac(ps_hell, weighted = TRUE, normalized = TRUE, 
                         parallel = FALSE, fast = TRUE)

# PERMANOVA test for weighted UniFrac ----
set.seed(4104)
permanova_wunifrac <- adonis2(dist_wunifrac ~ layer, data = metadata, permutations = 999)

# Test for homogeneity of dispersions ----
dispersion_wunifrac <- betadisper(dist_wunifrac, metadata$layer)
dispersion_test_wunifrac <- anova(dispersion_wunifrac)

# Convert distance matrix to dataframe ----
dist_wunifrac_matrix <- as.matrix(dist_wunifrac)
dist_wunifrac_df <- data.frame(
  Sample1 = rep(rownames(dist_wunifrac_matrix), each = ncol(dist_wunifrac_matrix)),
  Sample2 = rep(colnames(dist_wunifrac_matrix), times = nrow(dist_wunifrac_matrix)),
  WUniFrac = as.vector(dist_wunifrac_matrix)
)

# Remove self-comparisons
dist_wunifrac_df <- dist_wunifrac_df %>% filter(Sample1 != Sample2)

# Add comparison type columns ----
dist_wunifrac_df <- dist_wunifrac_df %>%
  mutate(
    Layer1 = extract_layer(Sample1),
    Layer2 = extract_layer(Sample2),
    ComparisonType = case_when(
      Layer1 == "Surface" & Layer2 == "Surface" ~ "Surface-Surface",
      Layer1 == "Bulk" & Layer2 == "Bulk" ~ "Bulk-Bulk",
      TRUE ~ "Mixed"
    ),
    Site1 = extract_site(Sample1),
    Site2 = extract_site(Sample2),
    ComparisonType2 = ifelse(Site1 == Site2, "Intra-site", "Inter-site")
  )

# Filter for inter-site comparisons only ----
dist_wunifrac_inter <- dist_wunifrac_df %>%
  filter(ComparisonType2 == "Inter-site") %>%
  filter(ComparisonType != "Mixed")

# Set factor levels
dist_wunifrac_inter$ComparisonType <- factor(dist_wunifrac_inter$ComparisonType, 
                                             levels = c("Surface-Surface", "Bulk-Bulk"))
# Statistical test: Wilcoxon test ----
wilcox_wunifrac <- dist_wunifrac_inter %>%
  pairwise_wilcox_test(WUniFrac ~ ComparisonType) %>%
  add_xy_position(x = "ComparisonType")
