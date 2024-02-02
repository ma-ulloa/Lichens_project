
#title: "Liquenes_análisis_16S"
#author: "Alejandra Ulloa"
#date: ' December 2023'

#This script calculates Alpha Diversity indeces for all lichen genus


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

#New otu-table
otu_rarefy <- otu.data.table(ps.rarefied)
otu_rarefy <- column_to_rownames(otu_rarefy, var = "Sample") 

##Metadata
meta_all <- sample.data.table(ps.rarefied)

#Estimate diversity indexes

data_richness <- estimateR(t(otu_rarefy))                                           # calculate richness and Chao1 using vegan package
data_evenness <- vegan::diversity(t(otu_rarefy)) / log(specnumber(t(otu_rarefy)))   # calculate evenness index using vegan package
data_shannon <- vegan::diversity(t(otu_rarefy), index = "shannon")                  # calculate Shannon index using vegan package
data_alphadiv <- cbind(meta_all, t(data_richness), data_shannon, data_evenness)     # combine all indices in one data table
write.table(data_alphadiv, file = "Meta-analysis/Alpha_diversity_table.txt", sep = "\t", dec = ",")
rm(data_richness, data_evenness, data_shannon)                                      # remove the unnecessary data/vector


## Stats Alpha-D

anova.host_genus = aov(data_shannon ~ host_genus, data = data_alphadiv)
summary(anova.host_genus)
TukeyHSD(anova.host_genus)

anova.host_genus.ev = aov(data_evenness ~ host_genus, data = data_alphadiv)
summary(anova.host_genus.ev)
TukeyHSD(anova.host_genus.ev)

anova.host_genus.rc = aov(S.ACE ~ host_genus, data = data_alphadiv)
summary(anova.host_genus.rc)
TukeyHSD(anova.host_genus.rc)

# Visualization

col_alp <- RColorBrewer::brewer.pal(3, "Dark2") 


P1 <- ggplot(data_alphadiv, aes(x=host_genus, y=S.chao1)) +
  geom_boxplot(fill=col_alp) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "A") +
  geom_point(size=0.6) +
  theme_classic() +
  geom_signif(comparisons = list(c("Stereocaulon", "Sticta"), c("Sticta", "Cora")),
              map_signif_level = TRUE,
              textsize = 4,
              margin_top = 0.08,
              step_increase = 0.08,
              tip_length = 0.01)

P2 <- ggplot(data_alphadiv, aes(x=host_genus, y=data_evenness)) +
  geom_boxplot(fill=col_alp) +
  labs(title= 'Evenness', x= ' ', y= '', tag = "B") +
  geom_point(size=0.6) + 
  theme_classic() +
  geom_signif(comparisons = list(c("Stereocaulon", "Cora"), c("Sticta", "Cora")),
              map_signif_level = TRUE,
              textsize = 4,
              margin_top = 0.08,
              step_increase = 0.08,
              tip_length = 0.01)

P3 <- ggplot(data_alphadiv, aes(x=host_genus, y=data_shannon)) +
  geom_boxplot(fill=col_alp) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "C") +
  geom_point(size=0.6) + 
  theme_classic() +
  geom_signif(comparisons = list(c("Stereocaulon", "Cora"), c("Sticta", "Cora")),
              map_signif_level = TRUE,
              textsize = 4,
              margin_top = 0.08,
              step_increase = 0.08,
              tip_length = 0.01)

# all plots together using the patchwork package


pdf("Meta-analysis/Alpha_diversity-meta-analysis.pdf",
    width = 9,
    height = 3.5)

(P1 | P2 | P3)

dev.off()