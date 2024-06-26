
#title: "Genus_Sticta"
#author: "Alejandra Ulloa"
#date: 'December 2023'

#This script was used to generate the relative abundance tables for different
#species of Sticta

#Load packages to be use

library(phyloseq)
library(data.table)
library(tidyr)
library(ggplot2)
library(ape)
library(vegan)
library(dplyr)
library(tibble)
library(phyloseqCompanion)
library(microViz)
library(microbiome)
library(ggsci)
library(purrr)
library(RColorBrewer)
library(patchwork)
library(ggsignif)
library(ggsignif)
library(tidyverse)

######## Import clean phyloseq object ##########

#load("Lichens_clean.RData")
load("Lichens_clean_review.RData")


# Keep only Sticta genus in Phyloseq object
physeq.sticta <-  subset_samples(ps.rarefied, host_genus == "Sticta")
physeq.sticta

#Fix tax table for Microviz

physeq.sticta  <- physeq.sticta  %>% tax_fix(min_length = 4,
                                             unknowns = c("g__Unknown_Family",
                                                          "f__Unknown_Family",
                                                          "g__Unknown_Family Genus",
                                                          "f__uncultured Family",
                                                          "g__uncultured",
                                                          "f__uncultured", 
                                                          "o__uncultured", 
                                                          "g__1174-901-12"),
                                             sep = " ", anon_unique = TRUE,
                                             suffix_rank = "classified")

#Prevalence 
lichen_melt = function(physeq,
                       includeSampleVars = character(),
                       omitZero = FALSE){
  # supports otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID] #Sum of counts
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. 
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}

mdt = lichen_melt(physeq.sticta)
lichen.prev = mdt[, list(Prevalence = mean(count > 0), #Mean of samples with counts greater than cero
                         TotalCounts = sum(count)),
                  by = TaxaID]
addPhylum = unique(copy(mdt[, list(TaxaID, Family)])) #Color taxa by Family
# Join by TaxaID
setkey(lichen.prev, TaxaID)
setkey(addPhylum, TaxaID)
lichen.prev <- addPhylum[lichen.prev]

# Detect Lichen core
core <- lichen.prev %>% filter(Prevalence >= 0.4) %>% pull(TaxaID)  # Take the ASVs that correspond to the core

# Select only the prevalent ASVs
physeq.sticta.core.P <- subset_taxa(physeq.sticta, ASVID %in% core)

# Generate an abundance table at different taxonomic levels
physeq.sticta.core <- physeq.sticta.core.P %>%
  aggregate_taxa(level = "Genus") %>%
  microbiome::transform(transform = "compositional") %>%
  psmelt() %>%
  filter(Genus != "") %>%
  select(Genus, Sample, Abundance)  %>%
  distinct() %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fn = list) %>%
  mutate(across(where(is.list), ~as.numeric(unlist(.)))) %>%  # Ensure values are numeric
  replace(is.na(.), 0)

# Ensure unique row names
physeq.sticta.core <- physeq.sticta.core %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()

# Convert to row names, transpose, and convert to data frame
physeq.sticta.core <- physeq.sticta.core %>%
  column_to_rownames(var = "Genus") %>%
  t() %>%
  as.data.frame()

physeq.sticta.core$sample.id <- rownames(physeq.sticta.core)

# Metadata
meta_st <- sample_data(physeq.sticta)  # Ensure to use sample_data() to get metadata as a data frame

# Merge OTU data with metadata
merged_data <- merge(meta_st, physeq.sticta.core, by = "sample.id")

# Convert microbial data to long format, excluding non-numeric columns
long_data <- merged_data %>%
  pivot_longer(cols = -c(sample.id, species.checked), names_to = "MicrobialSpecies", values_to = "frequency")

# Group by species.checked and MicrobialSpecies, and calculate median frequency
median_data <- long_data %>%
  group_by(species.checked, MicrobialSpecies) %>%
  summarise(median_frequency = median(frequency, na.rm = TRUE), .groups = 'drop')

# Filter out species with very few samples (less than 3 samples)
min_samples <- 3
species_sample_count <- merged_data %>%
  group_by(species.checked) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count >= min_samples)

median_data <- median_data %>%
  filter(species.checked %in% species_sample_count$species.checked)

# Convert data back to wide format
compiled_table <- median_data %>%
  pivot_wider(names_from = species.checked, values_from = median_frequency)

# Write the compiled table to a file
write.table(compiled_table, file = "Sticta/co-divergence/Reviewed/compiled_table_genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
