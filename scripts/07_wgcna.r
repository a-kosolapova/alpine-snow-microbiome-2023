# ==============================================================================
# Script 07: WGCNA - WEIGHTED GENE COEXPRESSION NETWORK ANALYSIS

# Load libraries ----
library(tidyverse)
library(phyloseq)
library(WGCNA)
library(DT)
library(htmlwidgets)

# Allow multi-threading
allowWGCNAThreads()

# Import filtered phyloseq object ----
ps <- readRDS("data/processed/phyloseq_filtered.rds")

# Hellinger transform ----
ps_hell <- transform(ps, transform = "hellinger")

# Extract OTU table and metadata ----
otu_hell <- t(as.matrix(otu_table(ps_hell)))
metadata <- data.frame(sample_data(ps_hell))

# Filter OTUs: keep only those with column sum > 0.05 ----
otu_filtered <- otu_hell[, colSums(otu_hell) > 0.05]

# Check for good samples and genes ----
gsg <- goodSamplesGenes(otu_filtered, verbose = 3)
if (!gsg$allOK) {
  cat("Removing problematic samples/genes\n")
  otu_filtered <- otu_filtered[gsg$goodSamples, gsg$goodGenes]
}

# Prepare trait data ----
traitData <- metadata %>%
  select(pH, fluoride, chloride, bromide, sulfate, phosphate, sodium, 
         potassium, magnesium, calcium, TIN, acids, elevation)

# Match samples between OTU and trait data ----
OTUSamples <- rownames(otu_filtered)
traitRows <- match(OTUSamples, rownames(traitData))
datTraits <- traitData[traitRows, ]

# Network construction ----------------------------------------------------

# Choose soft-thresholding power ----
powers <- c(1:10, seq(from = 11, to = 30, by = 1))
sft <- pickSoftThreshold(otu_filtered, powerVector = powers, verbose = 5, networkType = "signed")

# Use selected soft power (adjust based on your sft results)
softPower <- 9

# Calculate adjacency ----
adjacency <- adjacency(otu_filtered, power = softPower, type = "signed")

# Calculate TOM (Topological Overlap Matrix) ----
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Cluster taxa ----
TaxaTree <- hclust(as.dist(dissTOM), method = "average")

# Cut tree to identify modules ----
minModuleSize <- 100
dynamicMods <- cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)

# Calculate module eigengenes ----
MEList <- moduleEigengenes(otu_filtered, colors = dynamicColors)
MEs <- MEList$eigengenes

# Merge close modules ----
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
MEDissThres <- 0.90

merge <- mergeCloseModules(otu_filtered, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors <- merge$colors
MEs <- merge$newMEs

# Module-trait relationships ----------------------------------------------

nSamples <- nrow(otu_filtered)
MEs <- orderMEs(MEs)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create text matrix for heatmap ----
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

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
  "acids" = "'Org. Acids'",
  "TIN" = "'TIN'",
  "elevation" = "'Elevation'",
  "pH" = "pH"
)

convertedLabels <- sapply(colnames(datTraits), function(label) {
  if (label %in% names(label_conversion)) {
    parse(text = label_conversion[[label]])
  } else {
    label
  }
})

# Save module-trait heatmap ----
png(filename = "results/figures/wgcna_module_trait_heatmap.png", 
    width = 6000, height = 3500, res = 600)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = convertedLabels,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()


# Create module-taxonomy table --------------------------------------------

# Get OTU names from filtered data
otu_names <- colnames(otu_filtered)

# Create data frame with OTU and module assignments
otu_module_df <- data.frame(
  OTU = otu_names,
  moduleColor = moduleColors,
  stringsAsFactors = FALSE
)

# Get taxonomy from phyloseq
tax <- as.data.frame(tax_table(ps_hell))
tax_filtered <- tax[otu_names, , drop = FALSE]

# Combine module assignments with taxonomy
module_taxonomy <- cbind(otu_module_df, tax_filtered)

# Export module taxonomy table
write.csv(module_taxonomy, "results/wgcna_module_taxonomy.csv", row.names = FALSE)


# Enrichment analysis for selected modules --------------------------------

# Function to perform hypergeometric enrichment test for a module
perform_enrichment <- function(module_color, module_taxonomy, taxonomic_level = "Family") {
  
  # Rename taxonomy column
  df_annot <- module_taxonomy %>% rename(Organism = !!sym(taxonomic_level))
  
  # Total number of taxa
  Nt <- nrow(df_annot)
  
  # Count organisms in total dataset
  ModCountsM <- df_annot %>%
    count(Organism, name = "m")
  
  ModCountsN <- ModCountsM %>%
    mutate(n = Nt - m) %>%
    select(-m)
  
  # Filter for specific module
  df_module <- df_annot %>% filter(moduleColor == module_color)
  
  if (nrow(df_module) == 0) {
    return(NULL)
  }
  
  k <- nrow(df_module)  # Total taxa in module
  n0 <- length(unique(df_module$Organism))  # Unique organisms in module
  
  # Count organisms in module
  ModCountsQ <- df_module %>% count(Organism, name = "q")
  
  # Combine counts
  ModCounts <- left_join(ModCountsQ, ModCountsM, by = "Organism")
  ModCounts <- right_join(ModCountsN, ModCounts, by = "Organism")
  ModCounts <- ModCounts %>%
    mutate(k = k,
           p = phyper(q, m, n, k, lower.tail = FALSE),
           pAdj = p.adjust(p, n0, method = "BH"),
           Module = module_color)
  
  return(ModCounts)
}

# Select modules to analyze (adjust based on your module-trait heatmap)
selected_modules <- c("cyan", "darkgreen", "darkolivegreen", "darkred", 
                      "lightcyan", "sienna3", "skyblue3", "steelblue", "yellowgreen")

# Perform enrichment for all selected modules
enrichment_results <- list()
for (mod in selected_modules) {
  result <- perform_enrichment(mod, module_taxonomy, taxonomic_level = "Family")
  if (!is.null(result)) {
    enrichment_results[[mod]] <- result
  }
}

# Combine all enrichment results
ModFCounts <- bind_rows(enrichment_results)

# Filter for significant and abundant taxa
ModFCounts_sig <- ModFCounts %>% 
  filter(q > 9, pAdj < 0.05) %>%
  arrange(pAdj)

# Create interactive table
enrichment_table <- datatable(ModFCounts_sig, 
                              rownames = FALSE, 
                              filter = "top", 
                              options = list(order = list(list(3, 'asc'))))

# Save interactive table
saveWidget(enrichment_table, "results/wgcna_enrichment_family.html", selfcontained = TRUE)

# Export enrichment results as CSV
write.csv(ModFCounts_sig, "results/wgcna_enrichment_results.csv", row.names = FALSE)

# Save WGCNA workspace
save(moduleColors, MEs, moduleTraitCor, moduleTraitPvalue, module_taxonomy,
     file = "results/wgcna_workspace.RData")