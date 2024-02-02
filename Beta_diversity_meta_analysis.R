
#title: "Liquenes_análisis_16S"
#author: "Maria Alejandra Ulloa"

# This script calculates Beta-diversity metrics such as Bray-Curtis and unweighted
# Unifrac, produces respective ordination plots and creates a biplot with arrows
# Representing the taxa that displayed a greater association with the clustering.

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


#####

####################  Using Microviz for clustering of samples ##################

## Unifrac PCoA

ps.rarefied %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5) %>%
  ord_calc("PCoA") %>%
  ord_plot(color="host_genus", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Dark2") + theme_classic()

##Bray-Curtis PCoA

br.dist.b <- ps.rarefied %>%
  tax_transform("identity") %>%
  dist_calc(dist = "bray", binary = TRUE) %>%
  ord_calc("PCoA")

# Biplot 

#Taxonomy table
tax.rarefy <- as.data.frame(tax_table(ps.rarefied)) 

#Calculate distance
br_ord <- cmdscale(br.dist.b@dist, eig = TRUE, add = TRUE)
br.ordination <- as.data.frame(br_ord$points)

eig <- eigenvals(br_ord)
prop.x <- eig / sum(eig)
eig1 <- prop.x[1] %>% round(digits = 2) * 100
eig2 <- prop.x[2] %>% round(digits = 2) * 100

#Fit environmental variables
otu.fit <- otu_rarefy %>% t()

# Set seed for reproducibility
set.seed(123)

#Fit ASVs to ordination
fit.br <- envfit(br_ord, otu.fit, permutations = 999, na.rm = TRUE) #This might take a while... like 40min 

# Fit arrows to ordination space
en_coord_cat = as.data.frame(scores(fit.br, "vectors")) * ordiArrowMul(fit.br) 
en_coord_cat$ASVID <- rownames(en_coord_cat)

# Take 
arrows.br <- data.frame(R=fit.br$vectors$r, P=fit.br$vectors$pvals)
arrows.br$ASVID <- rownames(arrows.br)

# filter p < 0.01
arrows.p.br<- arrows.br[arrows.br$P<0.01,]

#From those select the ones with a magnitude > 0.28
arrows.p.br <- arrows.p.br[arrows.p.br$R>0.31,] 

arrows <- left_join(arrows.p.br, en_coord_cat, by= "ASVID")
arrows <- left_join(arrows, tax.rarefy, by="ASVID")

###Sample_data

data.sp <- sample.data.frame(ps.rarefied)

###plot dataframe
br_plot <- cbind(br.ordination, data.sp)

#Colors
colors <- c("#1B9E77", "#D95F02","#7570B3")

bray_curtis_plot <- br_plot %>% 
                    ggplot() + 
                    geom_point(aes(x=V1, y=V2,color = host_genus, shape = site)) +
                    scale_colour_manual(values = colors) +
                    xlab(paste("PCoA [", eig1, "%]")) +
                    ylab(paste("PCoA [", eig2, "%]")) +
                    theme_classic() 

br_biplot <- br_plot  %>% 
             ggplot() +
             geom_point(aes(x=V1, y=V2, color = host_genus, shape = site)) +
             scale_colour_manual(values = colors) +
             geom_segment(aes(x = 0,
                              y = 0,
                              xend = Dim1, 
                              yend = Dim2), 
                              data = arrows,
                              size =1, 
                              alpha = 0.8, 
                              colour = "grey30") +
             geom_text(aes(x = Dim1,
                           y = Dim2, 
                           label = Family), 
                           data = arrows, 
                           check_overlap = TRUE) +
             theme_minimal()

pdf("Meta-analysis/Bray-curtis-Biplot.pdf",
    width = 10,
    height = 5
)

br_biplot

dev.off()


################# In-group variance ############################

##Dispersion of samples

#Calculate
br.dist <- br.dist.b@dist
betadisper(br.dist, data.sp$host_genus)
permutest(dispersion.br)
anova(dispersion.br)
TukeyHSD(dispersion.br)
