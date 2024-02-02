
#title: "Genus_Sticta"
#author: "Alejandra Ulloa"
#date: 'December 2023'

#This script was used to analyze samples of the genus Sticta

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


################################################################################
########################### ALPHA DIVERSITY ####################################
################################################################################

##OTU TABLE
otu_st <- otu.data.table(physeq.sticta)
otu_st <- column_to_rownames(otu_st, var = "Sample") 

##Metadata

meta_st <- sample.data.table(physeq.sticta)

#Estimate diversity indexes

data_richness <- estimateR(t(otu_st))                                            # calculate richness and Chao1 using vegan package
data_evenness <- vegan::diversity(t(otu_st)) / log(specnumber(t(otu_st)))        # calculate evenness index using vegan package
data_shannon <- vegan::diversity(t(otu_st), index = "shannon")                   #Calculate Shannon index using vegan package
data_alphadiv <- cbind(meta_st, t(data_richness), data_shannon, data_evenness)   #Combine all indices in one data table
rm(data_richness, data_evenness, data_shannon)                                   #Remove the unnecessary data/vector
data_alphadiv <- data_alphadiv[-23,]

col_alp <- RColorBrewer::brewer.pal(4, "Dark2")

# Visualization

P1 <- ggplot(data_alphadiv, aes(x=host_species, y=S.chao1)) +
      geom_boxplot(fill=col_alp) +
      labs(title= 'Chao1',
           x= ' ', 
           y= '', 
           tag = "A") +
      geom_point() + 
      theme_classic() +
      geom_signif(comparisons = list(c("humboldtii", "andina"), 
                                     c("gyalocarpa","humboldtii")),
                  map_signif_level = TRUE,
                  textsize = 4,
                  margin_top = 0.08,
                  step_increase = 0.08,
                  tip_length = 0.01)

P2 <- ggplot(data_alphadiv, aes(x=host_species, y=data_shannon)) +
      geom_boxplot(fill=col_alp) +
      labs(title= 'Shannon',
           x= ' ', y= '',
           tag = "B") +
      geom_point() + 
      theme_classic() +
      geom_signif(comparisons = list(c("humboldtii", "andina"),
                                     c("gyalocarpa","humboldtii")),
                  map_signif_level = TRUE,
                  textsize = 4,
                  margin_top = 0.08,
                  step_increase = 0.08,
                  tip_length = 0.01)

# All plots together using the patchwork package

pdf("Sticta/Sticta_Alpha_D.pdf",
    width = 8 , 
    height = 5)

P1 | P2  + plot_annotation(tag_levels = 'A')

dev.off()

## Stats Alpha-D

anova.site = aov(data_shannon ~ site, data = data_alphadiv)
summary(anova.site)
TukeyHSD(anova.site)

anova.specie = aov(data_shannon ~ host_species, data = data_alphadiv)
summary(anova.specie)
TukeyHSD(anova.specie)

anova.specie.chao = aov(S.chao1 ~ host_species, data = data_alphadiv)
summary(anova.specie.chao)
TukeyHSD(anova.specie.chao)


########################################
############## Ordination ##############
########################################

# Unifrac PcoA

pdf("Sticta/PCoA-Sticta-Unifrac.pdf",
    width = 8 , 
    height = 5)

physeq.sticta %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc("PCoA") %>%
  ord_plot(color="host_species", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Dark2") + 
  stat_ellipse(aes(colour = host_species), linewidth = 0.2) +
  theme_classic()

dev.off()

# Bray-Curtis PCoA

pdf("Sticta/PCoA-Sticta-BC.pdf",
    width = 8 , 
    height = 5)

physeq.sticta %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc(dist = "bray", binary = TRUE) %>%
  ord_calc("PCoA") %>%
  ord_plot(color="host_species", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Dark2")  + 
  stat_ellipse(aes(colour = host_species), linewidth = 0.2) +
  theme_classic()

dev.off()

################ Statistical Analysis ###########

# Sample metadata
sample.st <- sample.data.frame(physeq.sticta)

br.dist.st <- physeq.sticta %>%
              tax_transform("identity", rank = "unique") %>%
              dist_calc(dist = "bray", binary = TRUE)
st.dist.br <- br.dist.st@dist

dispersion.br <- betadisper(st.dist.br, sample.st$host_species) %>%
                 permutest()


#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.083409 0.0278029 13.266    999  0.001 ***
#  Residuals 48 0.100597 0.0020958                         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

iav <- anova(dispersion.br)

#Analysis of Variance Table
#Response: Distances
#Df   Sum Sq   Mean Sq F value      Pr(>F)    
#Groups     3 0.083409 0.0278029  13.266 0.000001966 ***
#  Residuals 48 0.100597 0.0020958                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tukey.br <- TukeyHSD(dispersion.br)
#tukey.br 
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)
#$group
#diff          lwr          upr     p adj
#gyalocarpa-andina     -0.062770171 -0.116768576 -0.008771767 0.0167820
#humboldtii-andina     -0.099071719 -0.151828615 -0.046314823 0.0000470
#impressula-andina     -0.009824369 -0.063822773  0.044174036 0.9622419
#humboldtii-gyalocarpa -0.036301548 -0.080889263  0.008286167 0.1472739
#impressula-gyalocarpa  0.052945802  0.006895808  0.098995796 0.0183568
#impressula-humboldtii  0.089247350  0.044659635  0.133835065 0.0000153

########



#################################################################################
###### Rank of abundance  #######################

#Calculate the relative abundance
physeq_ra_st <- transform_sample_counts(physeq.sticta, function(x) x / sum(x))
physeq.ab.st <- psmelt(physeq_ra_st)
physeq.ab.st <- filter(physeq.ab.st, Abundance > 0)

#This is where the mean is calculated and the taxa to display is chosen
cluster_ps_st <- aggregate(Abundance ~ OTU + Family ,data=physeq.ab.st, mean)

# filtering and picking the number to display
cluster_ps_st = cluster_ps_st[order(-cluster_ps_st$Abundance),][1:1000,]

ggplot(cluster_ps_st, aes(x = reorder(OTU, -Abundance),y = Abundance)) +
  geom_point(aes(color = Family), size = 3) + 
  xlab("Rank") + 
  #scale_y_continuous(breaks = c(0.01:0.3)) +
  theme_bw()

#We are using all taxa > 0.01

ps.sticta.new <- physeq.ab.st %>% 
                 select(OTU, Sample, Abundance) %>% 
                 subset(Abundance >= 0.01) %>%
                 pivot_wider(names_from = Sample, values_from = Abundance)

##  [341 ASVs were selected from]

#Change NULL values to 0
ps.sticta.new[is.na(ps.sticta.new)] <- 0
ps.sticta.new <- column_to_rownames(ps.sticta.new, var = "OTU")
ps.sticta.new <- as.data.frame(ps.sticta.new)

