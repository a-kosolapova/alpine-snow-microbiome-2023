# ==============================================================================
# Script 06: DIFFERENTIAL ABUNDANCE ANALYSIS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Extract and prepare metadata ----
metadata <- data.frame(sample_data(ps))

# Log-transform chemical variables ----
log_transform <- function(x) {
  log(x + (min(x[x > 0]) / 2))
}

metadata$logTIN <- log_transform(metadata$TIN)
metadata$logCalcium <- log_transform(metadata$calcium)

# Update phyloseq object with transformed metadata ----
sample_data(ps) <- sample_data(metadata)

# Run ANCOM-BC2 with TIN and Calcium ----
set.seed(123)

output_tin_ca <- ancombc2(
  data = ps, 
  tax_level = "Genus",
  fix_formula = "logTIN + logCalcium", 
  rand_formula = NULL, 
  p_adj_method = "holm", 
  prv_cut = 0.1, 
  lib_cut = 1000, 
  s0_perc = 0.05
)

# Extract results ----
results <- output_tin_ca$res

# Process TIN results ----
df_tin <- results %>%
  select(taxon, lfc_logTIN, se_logTIN, diff_logTIN, p_logTIN) %>%
  filter(diff_logTIN == 1) %>%
  filter(grepl("^Genus:", taxon)) %>% 
  arrange(desc(lfc_logTIN)) %>%
  mutate(
    direct = ifelse(lfc_logTIN > 0, "Positive LFC", "Negative LFC"),
    fold_change = exp(lfc_logTIN),
    genus_name = gsub("Genus:|_[0-9]+", "", taxon),
    genus_name = gsub("_", " ", genus_name),
    lfc_lower = lfc_logTIN - 1.96 * se_logTIN,
    lfc_upper = lfc_logTIN + 1.96 * se_logTIN
  )

df_tin$taxon <- factor(df_tin$taxon, levels = df_tin$taxon)
df_tin$direct <- factor(df_tin$direct, levels = c("Positive LFC", "Negative LFC"))

# Process calcium results ----
df_calcium <- results %>%
  select(taxon, lfc_logCalcium, se_logCalcium, diff_logCalcium, p_logCalcium) %>%
  filter(diff_logCalcium == 1) %>%
  filter(grepl("^Genus:", taxon)) %>% 
  arrange(desc(lfc_logCalcium)) %>%
  mutate(
    direct = ifelse(lfc_logCalcium > 0, "Positive LFC", "Negative LFC"),
    fold_change = exp(lfc_logCalcium),
    genus_name = gsub("Genus:|_[0-9]+", "", taxon),
    genus_name = gsub("_", " ", genus_name),
    lfc_lower = lfc_logCalcium - 1.96 * se_logCalcium,
    lfc_upper = lfc_logCalcium + 1.96 * se_logCalcium
  )

df_calcium$taxon <- factor(df_calcium$taxon, levels = df_calcium$taxon)
df_calcium$direct <- factor(df_calcium$direct, levels = c("Positive LFC", "Negative LFC"))

# TIN plot ----
p_tin <- ggplot(df_tin, aes(x = reorder(genus_name, fold_change), 
                            y = fold_change)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = direct), size = 3) +
  geom_errorbar(aes(ymin = exp(lfc_lower),  
                    ymax = exp(lfc_upper),  
                    color = direct), 
                width = 0.2, alpha = 0.7) +
  coord_flip() +
  scale_y_log10() +
  scale_color_manual(values = c("Positive LFC" = "#d73027", 
                                "Negative LFC" = "#4575b4")) +
  labs(
    title = "Differential Abundance Response to TIN",
    x = "Genus",
    y = "Fold Change",
    color = "Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Calcium plot ----
p_calcium <- ggplot(df_calcium, aes(x = reorder(genus_name, fold_change), 
                                    y = fold_change)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = direct), size = 3) +
  geom_errorbar(aes(ymin = exp(lfc_lower),   
                    ymax = exp(lfc_upper),  
                    color = direct), 
                width = 0.2, alpha = 0.7) +
  coord_flip() +
  scale_y_log10() +
  scale_color_manual(values = c("Positive LFC" = "#d73027", 
                                "Negative LFC" = "#4575b4")) +
  labs(
    title = "Differential Abundance Response to Calcium",
    x = "Genus",
    y = "Fold Change", 
    color = "Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save plots ----
ggsave("results/figures/diff_abundance_tin.png", 
       plot = p_tin, 
       height = 4, width = 7, bg = "white", dpi = 600)

ggsave("results/figures/diff_abundance_calcium.png", 
       plot = p_calcium, 
       height = 3, width = 6, bg = "white", dpi = 600)

# Export results tables ----
write.csv(df_tin, "results/tables/diff_abundance_tin.csv", row.names = FALSE)
write.csv(df_calcium, "results/tables/diff_abundance_calcium.csv", row.names = FALSE)
