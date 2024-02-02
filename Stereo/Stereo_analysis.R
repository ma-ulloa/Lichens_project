
#title: "Genus_Stereocaulon"
#author: "Alejandra Ulloa"
#date: 'December 2023'

#This script was used to analyze samples of the genus Stereocualon 

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
library(geosphere)


# Import phyloseq object 
load("Lichens_clean.RData")

# Create a phyloseq subset for the genus Sterocaulon
physeq.ste <-  subset_samples(ps.rarefied, host_genus == "Stereocaulon")

########### ALPHA DIVERSITY #################

# ASV TABLE
otu_ste <- otu.data.table(physeq.ste)
otu_ste <- column_to_rownames(otu_ste, var = "Sample") 

# Metadata
meta_ste <- sample.data.table(physeq.ste)

# Estimate diversity indexes

data_richness <- estimateR(t(otu_ste))                                          # calculate richness and Chao1 using vegan package
data_evenness <- vegan::diversity(t(otu_ste)) / log(specnumber(t(otu_ste)))     # calculate evenness index using vegan package
data_shannon <- vegan::diversity(t(otu_ste), index = "shannon")                 # calculate Shannon index using vegan package
data_alphadiv <- cbind(meta_ste, t(data_richness), data_shannon, data_evenness) # combine all indices in one data table
rm(data_richness, data_evenness, data_shannon)  

# There is only one samples that was classified as Stereocaulon by morphology
data_alphadiv <- data_alphadiv[-13,]

# Visualization

col_alp <- RColorBrewer::brewer.pal(3, "Dark2")

comparison <- list(c("ramulosum", "pomiferum"),    
                   c("tomentosum", "pomiferum"),    
                   c("tomentosum", "ramulosum"))

P1 <- ggplot(data_alphadiv, aes(x=host_species, y=S.chao1)) +
      geom_boxplot() +
      labs(title= 'Chao1', x= ' ', y= '', tag = "A") +
      geom_point() + 
      theme_classic() +
      geom_signif(comparisons = comparison,
                  map_signif_level = TRUE,
                  textsize = 4,
                  margin_top = 0.08,
                  step_increase = 0.08,
                  tip_length = 0.01)

P2 <- ggplot(data_alphadiv, aes(x=host_species, y=data_shannon)) +
      geom_boxplot() +
      labs(title= 'Shannon', x= ' ', y= '', tag = "B") +
      geom_point() + 
      theme_classic() +
      geom_signif(comparisons = comparison,
                  map_signif_level = TRUE,
                  textsize = 4,
                  margin_top = 0.08,
                  step_increase = 0.08,
                  tip_length = 0.01)

pdf("Stereo/AlphaD_stereo.pdf",
    width = 14,
    height = 6
)

P1 | P2  + plot_annotation(tag_levels = 'A')

dev.off()

# Stats Alpha-D

# Shannon
anova.site = aov(data_shannon ~ host_species, data = data_alphadiv)
summary(anova.site)
TukeyHSD(anova.site)

# Chao1
anova.site.chao = aov(S.chao1 ~ host_species, data = data_alphadiv)
summary(anova.site)
TukeyHSD(anova.site.chao)

############################################

## Subset for Sterocaulon Atlanticus 

st.atlanticus <- c("Sample5", "Sample13", "Sample26",
                   "Sample27","Sample28", "Sample45",
                   "Sample51", "Sample80", "Sample85",
                   "Sample86","Sample89", "Sample89", 
                   "Sample90", "Sample91", "Sample92",
                   "Sample93","Sample95", "Sample96",
                   "Sample97", "Sample99", "Sample100", 
                   "Sample139", "Sample140")

ps.stereo <- subset_samples(physeq.ste, sample.id %in% st.atlanticus)


################################################################################
############################### Geospatial analysis ############################
################################################################################

#Longitude and latitude

sample_ste <- sample.data.frame(ps.stereo) 
cord_table <- data.frame(Sample.id= sample_ste$sample.id, 
                         Longitude=sample_ste$longitude,
                         Latitude= sample_ste$latitude) %>%
  column_to_rownames(var = "Sample.id") %>% 
  mutate_at(c('Longitude','Latitude'), as.numeric)

#Calculate the Haversine distance
dist.geo <- distm(cord_table, fun = distHaversine)
rownames(dist.geo) <- sample_ste$sample.id
colnames(dist.geo) <- sample_ste$sample.id
dist.geo <- as.matrix(dist.geo)

#Abundance matrix with bray-curtis
br.dist.stereo <- ps.stereo %>%
  dist_calc(dist = "bray", binary = TRUE)

df.dist.br <- br.dist.stereo@dist %>% as.matrix()

## Calculate mantel test for both matrices

#abundance vs geographic 
abund_geo  <- mantel(df.dist.br, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo

#library(MASS)
#write.matrix(dist.geo, file ="geospatial_dist_matrix.txt")

geo.long <- dist.geo %>%
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to = "B", values_to = "geo_dist")

#write.matrix(df.dist.br, file ="Bray_Curtis_matrix.txt")

br.long <- df.dist.br %>%
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to = "B", values_to = "bray")

## Plot correlation
sample_ste$A <- sample_ste$sample.id
sample_ste$B <- sample_ste$sample.id
df.geo.plot <- cbind(geo.long, bray.c = br.long$bray)
df.geo.plot$geo_dist <- df.geo.plot$geo_dist/1000  # Distance to Km
# Include collection site A 
df.geo.plot <- merge(df.geo.plot, 
                     sample_ste[, c("A", "site")],
                     by = "A", 
                     all.x = TRUE)
#Change rownames
colnames(df.geo.plot) <- c("A", "B", "geo_dist", "bray.c", "site_A")

# Include collection site B
df.geo.plot <- merge(df.geo.plot, 
                     sample_ste[, c("B", "site")],
                     by = "B", 
                     all.x = TRUE)
#Change rownames
colnames(df.geo.plot) <- c("A", "B", "geo_dist", "bray.c", "site_A", "site_B")

#Merge comparison sites
df.geo.plot$sites_comparison <- paste(df.geo.plot$site_A,"-",df.geo.plot$site_B)
df.geo.plot <- df.geo.plot[, !(names(df.geo.plot) %in% c("site_A", "site_B"))]
write.table(df.geo.plot, "Stereo/Geographical_dist_df.txt", sep = "\t", row.names = FALSE)

# Visualization
pdf("Stereo/Correlation_geo.pdf",
    width = 10,
    height = 5)

cmap = c("#7F3B08", "#33A02C", "#E08214", 
         "#E6F598", "#FEE0B6", "#D8DAEB", 
         "#B2ABD2", "#ABDDA4", "#542788", 
         "#2D004B", "#9E0142", "#D53E4F",
         "#C7EAE5", "#01665E", "#3288BD",
         "#66C2A5", "#5E4FA2" ) 

df.geo.plot %>%
  #filter(sample.id < B) %>%
  ggplot(aes(x=geo_dist, y=bray.c, color = sites_comparison)) +
  geom_point() +
  scale_color_manual(values = cmap) +
  theme_classic() +
  ylab("Bray-Curtis") + 
  xlab("Haversine distance (Km)") +
  ggtitle("Stereocaulon Atlanticum") +
  labs(color = "Collection sites comparison") +
  geom_smooth(method='lm',formula=y~x, se = TRUE) 

dev.off()

# Model for the association
model.stereo <- lm(bray.c ~ geo_dist, data = df.geo.plot)
summary <- summary(model.stereo)