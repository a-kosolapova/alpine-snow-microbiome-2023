# ==============================================================================
# Script 01: ANALYSIS OF THE CHEMICAL DATA

# Fig.2a: Paired differences in chemical content between layers ------------

# Load libraries ----
library(tidyverse)
library(reshape2)
library(ggplot2)

# Import filtered metadata ----
metadata <- read.csv("data/processed/metadata_filtered.csv")
rownames(metadata) <- metadata$Sample_ID
rows_to_remove <- c("H18_1_M", "H18_2_T_2")
metadata <- metadata[!(metadata$Sample_ID %in% rows_to_remove), ]

# Add subsite column ----
metadata$Subsite <- sub("_[^_]+$", "", metadata$Sample_ID)

# Log-transform chemical variables ----
# Function to handle zeros in log transformation
log_transform <- function(x) {
  log(x + (min(x[x > 0]) / 2))
}

# Apply log transformations
metadata$logChloride <- log_transform(metadata$chloride)
metadata$logTIN <- log_transform(metadata$TIN)
metadata$logAcids <- log_transform(metadata$acids)
metadata$logCalcium <- log_transform(metadata$calcium)
metadata$logSulfate <- log_transform(metadata$sulfate)
metadata$logFluoride <- log_transform(metadata$fluoride)
metadata$logBromide <- log_transform(metadata$bromide)
metadata$logPhosphate <- log_transform(metadata$phosphate)
metadata$logMagnesium <- log_transform(metadata$magnesium)
metadata$logSodium <- log_transform(metadata$sodium)
metadata$logPotassium <- log_transform(metadata$potassium)

# Calculate layer differences (Surface - Bulk) ----
env_alpha_data <- metadata %>% 
  select(Subsite, layer, logChloride, logTIN, logAcids, logCalcium, logSulfate,
         pH, logFluoride, logBromide, logPhosphate, logMagnesium, logSodium,
         logPotassium) %>% 
  group_by(Subsite, layer) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

env_delta <- env_alpha_data %>%
  pivot_wider(names_from = layer, 
              values_from = c(logChloride, logTIN, logAcids, logCalcium, logSulfate,
                              pH, logFluoride, logBromide, logPhosphate, logMagnesium, 
                              logSodium, logPotassium)) %>%
  mutate(across(ends_with("_Surface"), 
                ~ . - get(sub("_Surface", "_Bulk", cur_column())), 
                .names = "delta_{.col}")) %>% 
  select(Subsite, starts_with("delta_")) %>% 
  rename_with(~ gsub("_Surface$", "", .x)) %>%
  na.omit()

# Statistical testing: Wilcoxon test for each variable ----
variables <- c("logChloride", "logTIN", "logAcids", "logCalcium", "logSulfate",
               "pH", "logFluoride", "logBromide", "logPhosphate", "logMagnesium",
               "logSodium", "logPotassium")

pvals <- sapply(variables, function(var) {
  wilcox.test(env_delta[[paste0("delta_", var)]])$p.value
})

# Adjust p-values for multiple testing (FDR) ----
pvals_adj <- p.adjust(pvals, method = "fdr")

# Create significance labels ----
get_sig_label <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "n.s.")))
}

pvals_df <- data.frame(
  Variable = paste0("delta_", names(pvals)),
  p_adj = pvals_adj,
  sig_label = get_sig_label(pvals_adj)
)

# Prepare data for plotting ----
# Select variables to plot
plot_vars <- c("delta_pH", "delta_logTIN", "delta_logAcids", "delta_logPotassium", 
               "delta_logSodium", "delta_logCalcium", "delta_logMagnesium", 
               "delta_logChloride", "delta_logSulfate")

env_delta_long <- env_delta %>%
  select(Subsite, all_of(plot_vars)) %>%
  pivot_longer(cols = -Subsite, names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = plot_vars))

# Merge with p-values
env_delta_long <- merge(env_delta_long, pvals_df, by = "Variable")

# Define variable labels with chemical notation ----
labels <- c(
  "delta_pH" = "pH",
  "delta_logTIN" = "TIN",
  "delta_logAcids" = "Org.Acids",
  "delta_logPotassium" = expression(K^'+'),
  "delta_logSodium" = expression(Na^'+'),
  "delta_logCalcium" = expression(Ca^'2+'),
  "delta_logMagnesium" = expression(Mg^'2+'),
  "delta_logChloride" = expression(Cl^'-'),
  "delta_logSulfate" = expression(SO[4]^'2-')
)

# Create plot ----
fig1a <- env_delta_long %>%
  ggplot(aes(x = Value, y = Variable, fill = Value >= 0)) +
  geom_jitter(width = 0, height = 0.1, size = 2, shape = 21, color = "black", stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkblue", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "gray46", "FALSE" = "lightgray")) +
  theme_classic() +
  xlab(expression(Delta ~ "Layers (Surface - Bulk)")) +
  scale_y_discrete(limits = rev(names(labels)), labels = labels) +
  theme(
    legend.title = element_text(size = 14, face = "bold"), 
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 14),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25)
  ) + 
  geom_text(
    data = distinct(env_delta_long, Variable, sig_label),
    aes(x = Inf, y = Variable, label = sig_label),
    hjust = 1.1, size = 6, inherit.aes = FALSE
  )

# Save plot ----
ggsave("results/figures/fig2a.png", 
       plot = fig1a, 
       height = 4.5, 
       width = 6, 
       bg = "white", 
       dpi = 600)


# Fig. 2b: PCA  ----------------------------------------------------------

# Load libraries ----
library(tidyverse)
library(vegan)
library(ggrepel)

# Import filtered metadata ----
metadata <- read.csv("data/processed/metadata_filtered.csv")
rownames(metadata) <- metadata$Sample_ID

# Select chemical variables for PCA ----
chem_data <- metadata %>% 
  select(pH, sodium, potassium, magnesium, calcium, fluoride, chloride, 
         bromide, sulfate, phosphate, acids, TIN)

# Run PCA using rda with scaling ----
pca_chem <- rda(chem_data, scale = TRUE)

# Extract variance explained ----
eigenvalues <- pca_chem$CA$eig
variance_explained <- round(100 * eigenvalues / sum(eigenvalues), 2)

# Extract loadings (variable scores) ----
pca_variables <- as.data.frame(scores(pca_chem, display = "species"))
pca_variables$Variable <- rownames(pca_variables)

# Extract site scores ----
pca_scores <- as.data.frame(scores(pca_chem, display = "sites"))
pca_scores$layer <- metadata$layer
pca_scores$valley <- metadata$valley
pca_scores$site <- metadata$site

# Replace Other in valley with Peripheral
pca_scores$valley[pca_scores$valley == "Other"] <- "Peripheral"

# Calculate scaling factor for arrows ----
max_score <- max(abs(pca_scores$PC1), abs(pca_scores$PC2))
scale_factor <- max_score / max(abs(pca_variables$PC1), abs(pca_variables$PC2))

# Convert variable names to chemical notation ----
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
  "acids" = "'Organic\\nAcids'", 
  "TIN" = "'Total Inorganic Nitrogen'"
)

pca_variables <- pca_variables %>%
  mutate(Variable = label_conversion[Variable])

# Define color palette ----
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")

# Create PCA plot ----
fig2b <- ggplot(pca_scores, aes(x = PC1, y = PC2, shape = layer, color = valley)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(1, 19))+
  labs(x = "PC1", y = "PC2") +
  theme_minimal()  +
  geom_segment(data = pca_variables, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor), 
               arrow = arrow(length = unit(0.3, "cm")),
               color = "grey", size = 0.5)  +
  geom_text_repel(data = pca_variables,inherit.aes = FALSE, 
                  aes(x = ifelse(Variable == "phosphate", PC1 * scale_factor * 0.1, PC1 * scale_factor), 
                      y = ifelse(Variable == "phosphate", PC2 * scale_factor * 0.1, PC2 * scale_factor), label = Variable),
                  vjust = 1, hjust = 1, size = 4, color = "black", parse = TRUE) +
  labs(x = paste0("PC1 (", variance_explained[1], "%)"),
       y = paste0("PC2 (", variance_explained[2], "%)"),
       shape = "Snow\nLayer", 
       color = "Valley") +
  scale_color_manual(values = cbbPalette) 

# Save plot ----
ggsave("results/figures/fig2b.png", 
       plot = p, 
       height = 10, 
       width = 10, 
       bg = "white", 
       dpi = 600)


# Fig. 2c; Fig. 2d - Spearman correlation analysis for chemical parameters by layer -------------------------------------

# Load libraries ----
library(tidyverse)
library(corrplot)
library(Hmisc)
library(vegan)

# Import filtered metadata ----
metadata <- read.csv("data/processed/metadata_filtered.csv")
rownames(metadata) <- metadata$Sample_ID

# Select only variables used in correlation analysis ----
metadata_corr <- metadata %>% 
  select(layer, elevation, pH, sodium, potassium, magnesium, calcium, 
         fluoride, chloride, bromide, sulfate, phosphate, acids, TIN)

# Standardize numeric variables ----
metadata_corr_std <- cbind(
  layer = metadata_corr$layer,
  decostand(metadata_corr[, !names(metadata_corr) %in% "layer"], "standardize", 2)
)

# Separate by layer ----
# Surface layer
metadata_surface <- metadata_corr_std %>% 
  filter(layer == "Surface") %>%
  select(-layer)

# Bulk layer
metadata_bulk <- metadata_corr_std %>% 
  filter(layer == "Bulk") %>%
  select(-layer)

# Custom labels with chemical notation ----
custom_labels <- c(
  "sodium" = ":Na^'+'",
  "potassium" = ":K^'+'",
  "magnesium" = ":Mg^'2+'",
  "calcium" = ":Ca^'2+'",
  "fluoride" = ":F^'-'",
  "chloride" = ":Cl^'-'",
  "bromide" = ":Br^'-'",
  "sulfate" = ":SO[4]^'2-'",
  "phosphate" = ":PO[4]^'3-'",
  "acids" = ":'Organic Acids'", 
  "TIN" = "TIN",
  "pH" = "pH",
  "elevation" = "Elevation"
)

# Spearman correlation for Surface layer ----
cor_surface <- cor(metadata_surface, use = "complete.obs")
res_surface <- rcorr(as.matrix(metadata_surface), type = "spearman")
res_surface$Padj <- matrix(p.adjust(res_surface$P, method = "BH"), 
                           ncol = ncol(cor_surface))

# Apply custom labels
colnames(cor_surface) <- custom_labels[colnames(cor_surface)]
rownames(cor_surface) <- custom_labels[rownames(cor_surface)]

# Plot Surface layer correlation ----
png("results/figures/fig2c.png", 
    width = 2000, height = 2000, res = 300) 
corrplot(cor_surface, p.mat = res_surface$P, method = 'square', 
         diag = FALSE, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5, 
         tl.col = "black", insig = 'label_sig', pch.col = 'grey20')
dev.off()

# Spearman correlation for Bulk layer ----
cor_bulk <- cor(metadata_bulk, use = "complete.obs")
res_bulk <- rcorr(as.matrix(metadata_bulk), type = "spearman")
res_bulk$Padj <- matrix(p.adjust(res_bulk$P, method = "BH"), 
                        ncol = ncol(cor_bulk))

# Apply custom labels
colnames(cor_bulk) <- custom_labels[colnames(cor_bulk)]
rownames(cor_bulk) <- custom_labels[rownames(cor_bulk)]

# Plot Bulk layer correlation ----
png("results/figures/fig2d.png", 
    width = 2000, height = 2000, res = 300) 
corrplot(cor_bulk, p.mat = res_bulk$P, method = 'square', 
         diag = FALSE, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5, 
         tl.col = "black", insig = 'label_sig', pch.col = 'grey20')
dev.off()


