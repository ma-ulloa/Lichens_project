
#title: "Genus_Sticta"
#author: "Alejandra Ulloa"
#date: 'December 2023'

#This script was used to analyze samples of the genus Sticta with the updated
# taxonomic annotation.

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

#load("Lichens_clean.RData")
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
data_alphadiv <- data_alphadiv[-14,] # Ignore Cellulosa as it only has one sample

col_alp <- RColorBrewer::brewer.pal(8, "Spectral")

# Visualization

P1 <- ggplot(data_alphadiv, aes(x=species.checked, y=S.chao1)) +
      geom_boxplot(fill=col_alp) +
      labs(title= '',
           x= 'Specie', 
           y= 'Chao1', 
           tag = "A") +
      geom_point() +
      theme_classic()

P2 <- ggplot(data_alphadiv, aes(x=species.checked, y=data_shannon)) +
      geom_boxplot(fill=col_alp) +
      labs(title= '',
           x= 'Specie', 
           y= 'Shannon', 
           tag = "B") +
      geom_point() + 
      theme_classic()

# All plots together using the patchwork package

pdf("Sticta/Sticta_Alpha_D_reviewed.pdf",
    width = 12, 
    height = 10)

P1 / P2  + plot_annotation(tag_levels = 'A')

dev.off()

## Stats Alpha-D

anova.site = aov(data_shannon ~ site, data = data_alphadiv)
summary(anova.site)
TukeyHSD(anova.site)

anova.specie = aov(data_shannon ~ species.checked, data = data_alphadiv)
summary(anova.specie)
TukeyHSD(anova.specie)

anova.specie.chao = aov(S.chao1 ~ species.checked, data = data_alphadiv)
summary(anova.specie.chao)
TukeyHSD(anova.specie.chao)


########################################
############## Ordination ##############
########################################

# Unifrac PcoA

pdf("Sticta/PCoA-Sticta-Unifrac-review.pdf",
    width = 8 , 
    height = 5)

physeq.sticta %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc("PCoA") %>%
  ord_plot(color="species.checked", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Spectral") + 
  stat_ellipse(aes(colour = species.checked), linewidth = 0.2) +
  theme_classic()

dev.off()

# Bray-Curtis PCoA

pdf("Sticta/PCoA-Sticta-BC-review.pdf",
    width = 8 , 
    height = 5)

physeq.sticta %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc(dist = "bray", binary = TRUE) %>%
  ord_calc("PCoA") %>%
  ord_plot(color="species.checked", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Spectral")  + 
  stat_ellipse(aes(colour = species.checked), linewidth = 0.2) +
  theme_classic()

dev.off()

################ Statistical Analysis ###########

# Sample metadata
sample.st <- sample.data.frame(physeq.sticta)

br.dist.st <- physeq.sticta %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc(dist = "bray", binary = TRUE)
st.dist.br <- br.dist.st@dist

dispersion.br <- betadisper(st.dist.br, sample.st$species.checked)
dispersion.br.perm <- dispersion.br %>% permutest()

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     8 0.34326 0.042908 14.738    999  0.001 ***
#  Residuals 43 0.12519 0.002911                         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#iav <- anova(dispersion.br)

iav <- anova(dispersion.br)

#Analysis of Variance Table

#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     8 0.34326 0.042908  14.738 4.254e-10 ***
#  Residuals 43 0.12519 0.002911                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey.br <- TukeyHSD(dispersion.br)
# tukey.br
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#diff         lwr         upr     p adj
#Sticta aff gyalocarpa-Sticta aff cellulosa        0.4825332933  0.30020544  0.66486115 0.0000000
#Sticta andina-Sticta aff cellulosa                0.5400233611  0.35319290  0.72685382 0.0000000
#Sticta arachnofuliginosa-Sticta aff cellulosa     0.3167330677  0.10099984  0.53246630 0.0006199
#Sticta brevior-Sticta aff cellulosa               0.4827648425  0.29250602  0.67302367 0.0000000
#Sticta parahumboldtii-Sticta aff cellulosa        0.4109207129  0.22066189  0.60117954 0.0000004
#Sticta pseudohumboldtii-Sticta aff cellulosa      0.4205179031  0.23368745  0.60734836 0.0000001
#Sticta rhizinata-Sticta aff cellulosa             0.4613801912  0.25798495  0.66477543 0.0000001
#Sticta sylvatica-Sticta aff cellulosa             0.4860058273  0.28906923  0.68294242 0.0000000
#Sticta andina-Sticta aff gyalocarpa               0.0574900678 -0.02057806  0.13555819 0.3078638
#Sticta arachnofuliginosa-Sticta aff gyalocarpa   -0.1658002255 -0.29895367 -0.03264679 0.0056867
#Sticta brevior-Sticta aff gyalocarpa              0.0002315492 -0.08571863  0.08618173 1.0000000
#Sticta parahumboldtii-Sticta aff gyalocarpa      -0.0716125803 -0.15756276  0.01433760 0.1711836
#Sticta pseudohumboldtii-Sticta aff gyalocarpa    -0.0620153902 -0.14008351  0.01605273 0.2193593
#Sticta rhizinata-Sticta aff gyalocarpa           -0.0211531021 -0.13321842  0.09091222 0.9994172
#Sticta sylvatica-Sticta aff gyalocarpa            0.0034725341 -0.09639255  0.10333761 1.0000000
#Sticta arachnofuliginosa-Sticta andina           -0.2232902933 -0.36254549 -0.08403509 0.0001511
#Sticta brevior-Sticta andina                     -0.0572585186 -0.15238793  0.03787089 0.5747736
#Sticta parahumboldtii-Sticta andina              -0.1291026481 -0.22423206 -0.03397324 0.0019089
#Sticta pseudohumboldtii-Sticta andina            -0.1195054580 -0.20757818 -0.03143274 0.0019133
#Sticta rhizinata-Sticta andina                   -0.0786431699 -0.19789420  0.04060786 0.4530733
#Sticta sylvatica-Sticta andina                   -0.0540175337 -0.16188415  0.05384908 0.7805440
#Sticta brevior-Sticta arachnofuliginosa           0.1660317747  0.02220962  0.30985393 0.0132360
#Sticta parahumboldtii-Sticta arachnofuliginosa    0.0941876452 -0.04963451  0.23800980 0.4625055
#Sticta pseudohumboldtii-Sticta arachnofuliginosa  0.1037848354 -0.03547037  0.24304004 0.2930408
#Sticta rhizinata-Sticta arachnofuliginosa         0.1446471234 -0.01615093  0.30544518 0.1081009
#Sticta sylvatica-Sticta arachnofuliginosa         0.1692727596  0.01672633  0.32181919 0.0197621
#Sticta parahumboldtii-Sticta brevior             -0.0718441295 -0.17354175  0.02985349 0.3610767
#Sticta pseudohumboldtii-Sticta brevior           -0.0622469394 -0.15737635  0.03288247 0.4636413
#Sticta rhizinata-Sticta brevior                  -0.0213846513 -0.14593829  0.10316899 0.9997098
#Sticta sylvatica-Sticta brevior                   0.0032409849 -0.11046041  0.11694238 1.0000000
#Sticta pseudohumboldtii-Sticta parahumboldtii     0.0095971902 -0.08553222  0.10472660 0.9999950
#Sticta rhizinata-Sticta parahumboldtii            0.0504594782 -0.07409416  0.17501312 0.9189030
#Sticta sylvatica-Sticta parahumboldtii            0.0750851144 -0.03861628  0.18878651 0.4512397
#Sticta rhizinata-Sticta pseudohumboldtii          0.0408622881 -0.07838874  0.16011332 0.9680987
#Sticta sylvatica-Sticta pseudohumboldtii          0.0654879242 -0.04237869  0.17335454 0.5635711
#Sticta sylvatica-Sticta rhizinata                 0.0246256362 -0.10990767  0.15915894 0.9995345




