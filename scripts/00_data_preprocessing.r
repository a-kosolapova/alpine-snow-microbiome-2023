# ==============================================================================
# Script 00: DATA IMPORT AND FILTERING

# Load libraries ----
library(tidyverse)    
library(phyloseq)     
library(ape)          
library(decontam)     

# Helper functions ----
clean_taxonomy_prefixes <- function(tax_matrix) {
  gsub("d__", "",
       gsub("p__", "",
            gsub("o__", "",
                 gsub("c__", "",
                      gsub("g__", "",
                           gsub("s__", "",
                                gsub("f__", "", tax_matrix)))))))
}

remove_unwanted_taxa <- function(ps_obj) {
  ps_obj %>%
    subset_taxa((Kingdom != "Unassigned") | is.na(Kingdom)) %>%
    subset_taxa((Class != "Chloroplast") | is.na(Class)) %>%
    subset_taxa((Family != "Mitochondria") | is.na(Family)) %>%
    subset_taxa((Kingdom != "Archaea") | is.na(Kingdom))
}

# Import raw data ----
otu_raw <- read_tsv("data/raw/otu_table.tsv")
taxonomy_raw <- read.csv("data/raw/taxonomy_table.tsv")
metadata <- read.csv("data/raw/metadata_table.csv", row.names = 1)

# Filter taxonomy by confidence >= 0.7 ----
taxonomy_filtered <- taxonomy_raw %>% filter(Confidence >= 0.7)
otu_filtered <- otu_raw[otu_raw$`#OTU ID` %in% taxonomy_filtered$Cluster, ]
otu_matrix <- otu_filtered %>% 
  column_to_rownames(var = "#OTU ID") %>%
  as.matrix()

# Process taxonomy table ----
taxonomy_processed <- taxonomy_filtered %>% 
  column_to_rownames(var = "Cluster") %>%
  select(Taxonomy) %>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ") %>%
  as.matrix()

taxonomy_clean <- clean_taxonomy_prefixes(taxonomy_processed)

# Clean metadata: replace NA with 0 for sample rows ----
metadata_clean <- metadata %>%
  mutate(across(everything(), ~ ifelse(type == "sample" & is.na(.), 0, .)))

# Create initial phyloseq object ----
ps_initial <- phyloseq(
  otu_table(otu_matrix, taxa_are_rows = TRUE),
  tax_table(taxonomy_clean),
  sample_data(metadata_clean)
)

# Taxonomic filtering: remove unwanted taxa ----
ps_bacteria_only <- remove_unwanted_taxa(ps_initial)

# Decontamination ----
# Identify negative controls
sample_data(ps_bacteria_only)$is.neg <- sample_data(ps_bacteria_only)$type != "sample"

# Find contaminants using prevalence method
contaminants <- isContaminant(ps_bacteria_only, method = "prevalence", neg = "is.neg", threshold = 0.5)

# Remove identified contaminants
ps_decontaminated <- prune_taxa(!contaminants$contaminant, ps_bacteria_only)

# Manual removal of specific problematic taxa
ps_manual_cleaned <- prune_taxa(!taxa_names(ps_decontaminated) %in% "NC9_10.s1", ps_decontaminated)

# Remove control samples, keep only biological samples ----
ps_samples_only <- subset_samples(ps_manual_cleaned, type == "sample")

# Abundance-based filtering ----
# Remove singletons (total abundance = 1)
ps_no_singletons <- filter_taxa(ps_samples_only, function(x) sum(x) > 1, TRUE)

# Remove low-abundance taxa (total abundance < 10)
ps_min_abundance <- filter_taxa(ps_no_singletons, function(x) sum(x) >= 10, TRUE)

# Prevalence filtering: keep taxa present in >= 2 samples ----
otu_prevalence <- rowSums(otu_table(ps_min_abundance) > 0)
taxa_to_keep <- names(otu_prevalence[otu_prevalence >= 2])
ps_filtered <- prune_taxa(taxa_to_keep, ps_min_abundance)

# Add phylogenetic tree ----
phylo_tree <- read.tree("data/processed/rooted_tree.nwk")
ps_final <- merge_phyloseq(ps_filtered, phylo_tree)

# Export filtered data ----
# Create output directory if it doesn't exist ----
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

# Export OTU table ----
otu_table_final <- as.data.frame(otu_table(ps_final))
otu_table_final$OTU_ID <- rownames(otu_table_final)
otu_table_final <- otu_table_final[, c("OTU_ID", setdiff(names(otu_table_final), "OTU_ID"))]
write.csv(otu_table_final, "data/processed/otu_table_filtered.csv", row.names = FALSE)

# Export taxonomy table  ----
taxonomy_table_final <- as.data.frame(tax_table(ps_final))
taxonomy_table_final$OTU_ID <- rownames(taxonomy_table_final)
taxonomy_table_final <- taxonomy_table_final[, c("OTU_ID", setdiff(names(taxonomy_table_final), "OTU_ID"))]
write.csv(taxonomy_table_final, "data/processed/taxonomy_table_filtered.csv", row.names = FALSE)

# Export metadata ----
metadata_final <- as.data.frame(sample_data(ps_final))
metadata_final$Sample_ID <- rownames(metadata_final)
metadata_final <- metadata_final[, c("Sample_ID", setdiff(names(metadata_final), "Sample_ID"))]
write.csv(metadata_final, "data/processed/metadata_filtered.csv", row.names = FALSE)

# Save final filtered phyloseq object ----
saveRDS(ps_final, "data/processed/phyloseq_filtered.rds")