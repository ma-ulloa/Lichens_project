
#This script builts the visualization for the core microbiome.
# It uses the lichen_melt() function from the script Prevalence. Available :
# https://github.com/mariaasierra/Lichen_Microbiome

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
library(ggsci)
library(purrr)
library(RColorBrewer)
library(patchwork)
library(ggsignif)

######## Import clean phyloseq object ##########

load("Lichens_clean.RData")

#######

######################## CORE PREVALENCE ######################

# Fix taxonomy labels for microviz to work
ps.asv <- ps.rarefied %>%
  tax_fix(min_length = 4,
          unknowns = c("g__Unknown_Family",
                       "f__Unknown_Family",
                       "g__Unknown_Family Genus",
                       "g__uncultured",
                       "f__uncultured", 
                       "o__uncultured", 
                       "g__1174-901-12"),
          sep = " ", anon_unique = TRUE,
          suffix_rank = "classified")

# Collapse to Family level
ps.family <- tax_agg(ps.asv, rank = "Family") 


##Function
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

mdt = lichen_melt(ps.family) # Change for ps.asv for asv
lichen.prev = mdt[, list(Prevalence = mean(count > 0), #Mean of samples with counts greater than cero
                         TotalCounts = sum(count)),
                  by = TaxaID]
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)])) #Color taxa by phylum
# Join by TaxaID
setkey(lichen.prev, TaxaID)
setkey(addPhylum, TaxaID)
lichen.prev <- addPhylum[lichen.prev]
showPhyla = lichen.prev[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(lichen.prev, Phylum)

#Color
color.prev <- ggsci::pal_futurama(palette = c("planetexpress"), alpha = 1)

pdf("Meta-analysis/Microbiome_core_collapsed.pdf",
    width = 10,
    height = 5
)

#Prevalence plot of total OTUs
ggplot(lichen.prev[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_jitter(size = 5, alpha = 0.7) + 
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  annotate(geom='label', x=0.5, y=9500, label="Pan", fontface=2) +
  annotate(geom='label', x=0.05, y=9500, label="Peripheral", fontface=2) +
  annotate(geom='label', x=1, y=9500, label="Core", fontface=2) +
  geom_vline(xintercept=0.25, color='gray20') +
  geom_vline(xintercept=0.9, color='gray20') +
  scale_y_log10() +  theme_minimal() + ylab("Total Counts") +
  theme(strip.text = element_text(size=20),
        axis.text.x = element_text(colour = "black", size=11 ), 
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.x = element_text(size=13), 
        axis.title.y = element_text(size=13),
        legend.text  = element_text(size=12, colour = "gray40"),
        legend.title = element_text(size = 14))

dev.off()
