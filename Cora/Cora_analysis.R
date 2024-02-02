
#title: "Genus_Sticta"
#author: "Alejandra Ulloa"
#date: 'December 2023'

#This script was used to analyze samples of the genus Cora

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
library(ggsignif)
library(tidyverse)

######## Import clean phyloseq object ##########

load("Lichens_clean.RData")


# Keep only Stereocaulon genus in Phyloseq object
physeq.cora <-  subset_samples(ps.rarefied, host_genus == "Cora")

########### ALPHA DIVERSITY #################

# OTU TABLE
otu.cora <- otu.data.table(physeq.cora)
otu.cora <- column_to_rownames(otu.cora, var = "Sample") 

# Metadata
meta_cora <- sample.data.table(physeq.cora)

# Estimate diversity indexes

data_richness <- estimateR(t(otu.cora))                                         # calculate richness and Chao1 using vegan package
data_evenness <- vegan::diversity(t(otu.cora)) / log(specnumber(t(otu.cora)))   # calculate evenness index using vegan package
data_shannon <- vegan::diversity(t(otu.cora), index = "shannon")                # calculate Shannon index using vegan package
data_alphadiv <- cbind(meta_cora, t(data_richness), data_shannon, data_evenness)# combine all indices in one data table
rm(data_richness, data_evenness, data_shannon)                                  # remove the unnecessary data/vector


#Ditch species with only one sample
data_alphadiv <- data_alphadiv[-c(17,18,30),]


# Visualization
col_alp <- RColorBrewer::brewer.pal(2, "Dark2") 

P1 <- ggplot(data_alphadiv, aes(x=host_species_ITS, y=S.chao1)) +
      geom_boxplot() +
      labs(title= 'Chao1', x= ' ', y= '', tag = "A") +
      geom_signif(comparisons = list(c("Cora_paraciferrii", "Cora_sp")),
              map_signif_level = TRUE,
              textsize = 4,
              margin_top = 0.08,
              step_increase = 0.08,
              tip_length = 0.01) +
      geom_point() + 
      theme_classic()

P2 <- ggplot(data_alphadiv, aes(x=host_species_ITS, y=data_shannon)) +
  geom_boxplot() +
  labs(title= 'Shannon', x= ' ', y= '', tag = "B") +
  geom_signif(comparisons = list(c("Cora_paraciferrii", "Cora_sp")),
              map_signif_level = TRUE,
              textsize = 4,
              margin_top = 0.08,
              step_increase = 0.08,
              tip_length = 0.01) +
  geom_point() + 
  theme_classic()

# all plots together using the patchwork package

pdf("Cora/alphaD_cora.pdf",
    width = 12,
    height = 6
)

P1 | P2  + plot_annotation(tag_levels = 'A')

dev.off()

## Stats Alpha-D

anova.site = aov(data_shannon ~ host_species_ITS, data = data_alphadiv)
TukeyHSD(anova.site)

anova.site.chao = aov(S.chao1 ~ host_species_ITS, data = data_alphadiv)
TukeyHSD(anova.site.chao)