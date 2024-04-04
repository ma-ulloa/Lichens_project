
#title: "Liquenes_análisis_16S"
#author: "Alejandra Ulloa"
#date: '2023-04-18'

## This script was used to make the feature selection of samples and ASVs
## to be used in the rest of the analysis. The idea was to remove NA, NULL
#values, ASVs corresponding to Archaea, Eukaryota and Unassigned


#Load packages to be used

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

######## Import objects for 16S data sets ##########

# ASV table
table_16S <- read.delim2("Raw_data/table-16S.tsv")
asv_counts <- column_to_rownames(table_16S, var = colnames(table_16S)[1])

# Taxonomy
taxmat_16S <- read.delim2("Raw_data/taxonomy-16S.tsv") 
taxmat_16S <- taxmat_16S %>% column_to_rownames(var = "Feature.ID")

# Sample metadata
sampledata_16S <- read.delim2("Raw_data/sample-metadata-16S-edit-stereo.tsv")
rownames(sampledata_16S) <- sampledata_16S$sample.id

# Phylogenetic tree - rooted
tree_16S = read.tree("Raw_data/tree-16S.nwk")


#################### Phyloseq ############################

# Create Phyloseq objects
OTU_16S = otu_table(asv_counts, taxa_are_rows = TRUE)
TAX_16S = tax_table(as.matrix(taxmat_16S))
Samples_16S = sample_data(sampledata_16S)

# Merge into phyloseq objects
physeq16s = phyloseq(OTU_16S, TAX_16S,Samples_16S)
# Add tree to phyloseq object
ps_16S = merge_phyloseq(physeq16s, tree_16S)


# Feature selection

unique(tax_table(ps_16S)[, "Domain"])
table(tax_table(ps_16S)[, "Domain"], exclude = NULL)

## Remove Archaea, Eukaryota and Unassigned

ps <- subset_taxa(ps_16S, !is.na(Domain) & !Domain %in%
                    c("Unassigned", "d__Eukaryota", "d__Archaea"))
      table(tax_table(ps)[, "Domain"], exclude = NULL)

## Filtering NA taxa or low abundant taxa phylum
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))

### Class
ps <- subset_taxa(ps, !is.na(Class) & !Class %in% c(""))

### Order
ps <- subset_taxa(ps, !is.na(Order) & !Order %in% c("","o__Chloroplast"))

###Family
ps <- subset_taxa(ps, !is.na(Family) & !Family %in% c("","f__Mitochondria"))

###Genus
ps <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c(""))

###### Update Phyloseq object with samples to remove ######

#First, create a list of the samples that you want to remove     
Samples_toRemove <- c("Sample11", "Sample16", "Sample25",
                      "Sample43", "Sample44", "Sample46",
                      "Sample50", "Sample52", "Sample56",
                      "Sample63", "Sample78", "Sample79", "Sample62",
                      "Sample81", "Sample87", "Sample88", "Sample32",
                      "Sample131", "Sample136", "Sample146", "Sample152")


# To see what samples get removed, run the following; note, I have a column called "SampleID"
subset_samples(ps, sample.id %in% Samples_toRemove)
# This will return a ps object that contains the samples you want to remove

# To remove those from your phyloseq object
ps <- subset_samples(ps, !(sample.id %in% Samples_toRemove))
# This will return a ps object with the samples removed

# Fix tax table for Microviz
ps <- ps %>% tax_fix()

# Add column to tax table
tax_table(ps) <- cbind(tax_table(ps),rownames(tax_table(ps)))

# Add ASVs to tax table as column
colnames(tax_table(ps)) <- 
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASVID")


# Rarefraction curves

sum_before <- sample_sums(ps)

#
otu_all <- otu.data.table(ps) %>% column_to_rownames(var = "Sample")

## Rarefraction curve
min_value <- min(rowSums(t(otu_all)))

## Histogram of read distribution

df <- data.frame(ASVs=rowSums(t(otu_table(ps))>0), reads=sample_sums(ps), sample_data(ps))

#write.table(df, "Before-rarefraction-read-data.txt", sep = ",")

hist.reads <- ggplot(df, aes(x=reads)) + 
  geom_histogram(bins=50, color='black', fill='grey') + 
  theme_bw()  +
  geom_vline(xintercept=1800, color= "red", linetype='dashed') +
  labs(title="Histogram: Reads per sample") +
  xlab("Read Count") + 
  ylab("Sample Count")

## Plot

rarecurve_data <- rarecurve(t(otu_all), 
                            step=100 , 
                            lwd=2, 
                            ylab="ASVs", 
                            label=F,
                            xlim=c(0, 5000),
                            ylim=c(0,300),
                            main="Rarefaction Curve for all samples")

plot.colors <- c(pal_aaas(palette= c("default"), alpha = 1)(10), 
                 pal_futurama(palette = c("planetexpress"), alpha = 1)(10),
                 pal_cosmic(palette = c("hallmarks_light"))(10),
                 pal_npg(palette = c("nrc"), alpha = 1)(10),
                 brewer.pal(11,"BrBG"),
                 brewer.pal(11, "PiYG"),
                 brewer.pal(11, "PRGn"),
                 brewer.pal(11, "PuOr"),
                 brewer.pal(11, "RdBu"),
                 brewer.pal(11, "RdGy"),
                 brewer.pal(11, "RdYlBu"),
                 brewer.pal(11, "RdYlGn")
)

pdf("Rarefraction-curves.pdf",
    width = 7,
    height = 5
)

map_dfr(rarecurve_data, bind_rows) %>%  
  bind_cols(Group = row.names(t(otu_all))) %>%
  pivot_longer(-Group) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  ggplot(aes(x=n_seqs, y=value, color =Group)) +
  geom_line(show.legend = FALSE) +
  scale_color_manual(values = plot.colors) +
  ggtitle("Rarefraction curves") +
  labs(x = "Sample size", y= "ASVs Count") +
  theme_bw()


dev.off()

# Rarefy

set.seed(1)  #For a reproducible code
ps.rarefied <- rarefy_even_depth(ps, rngseed=1, sample.size=1900, replace=F) 
#Sample117Sample138Sample54 were remove. 
#1061 OTUs were removed because they are no longer present in any sample after random subsampling


#New otu-table
otu_rarefy <- otu.data.table(ps.rarefied)
otu_rarefy <- column_to_rownames(otu_rarefy, var = "Sample") 

##Metadata
meta_all <- sample.data.table(ps.rarefied)

#Save data for next analysis
save(ps.rarefied,otu_rarefy, meta_all, file="Lichens_clean.RData")

