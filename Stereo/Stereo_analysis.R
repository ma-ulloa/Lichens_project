
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

sample_ste <- ps.stereo %>% sample.data.frame() %>%
              arrange(site)

#Join this to phyloseq object


cord_table <- data.frame(Sample.id= sample_ste$sample.id, 
                         Longitude=sample_ste$longitude,
                         Latitude= sample_ste$latitude) %>%
  column_to_rownames(var = "Sample.id") %>% 
  mutate_at(c('Longitude','Latitude'), as.numeric)

#Calculate the Haversine distance
dist.geo <- distm(cord_table, fun = distHaversine)
rownames(dist.geo) <- sample_ste$sample.id
colnames(dist.geo) <- sample_ste$sample.id
#dist.geo <- dist.geo/1000
dist.geo <- as.matrix(dist.geo)

#Abundance matrix with bray-curtis

br.dist.stereo <- ps.stereo %>% dist_calc(dist = "bray", binary = TRUE)

df.dist.br <- br.dist.stereo@dist %>% as.matrix()
df.dist.br <- df.dist.br[sample_ste$sample.id, sample_ste$sample.id]

## Calculate mantel test for both matrices

#abundance vs geographic 
abund_geo  <- mantel(df.dist.br, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo


## Correlogram of matrices

library(ComplexHeatmap)
library(circlize)
library(pheatmap)

# Generate annotations
anot_bottom <- HeatmapAnnotation(Collection_site = sample_ste$site, 
                                 col = list(Collection_site = c("Paramo de Chingaza" ="#1B9E77",
                                                                "Paramo de Sumapaz" = "#E6F598",
                                                                "Paramo del Verjon" = "#CA9822",
                                                                "Volcan de Cerro Machin" = "#351B9E")))
anot_right <- rowAnnotation(Collection_site = sample_ste$site, 
                            show_annotation_name = FALSE,
                            show_legend = FALSE,
                            col = list(Collection_site = c("Paramo de Chingaza" ="#1B9E77",
                                                           "Paramo de Sumapaz" = "#E6F598",
                                                           "Paramo del Verjon" = "#CA9822",
                                                           "Volcan de Cerro Machin" = "#351B9E")))

### Correlation plot

cor_geo_bray <- cor(df.dist.br,dist.geo)

# Correlation plot
nm <- rownames(cor_geo_bray)

# Color pallete
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#ffa500", "white", "#00a06d"))


pdf("Stereo/Corplot_geo_bray.pdf",
    width = 12,
    height = 8)

Heatmap(cor_geo_bray, name = "Correlation", 
        col = col_fun, 
        rect_gp = gpar(type = "none"), 
        bottom_annotation = anot_bottom,
        right_annotation = anot_right,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "grey", fill = NA))
          if(i > j) {
            grid.circle(x = x, y = y, r = abs(cor_geo_bray[i, j])/2 * min(unit.c(width, height)), 
                        gp = gpar(fill = col_fun(cor_geo_bray[i, j]), col = NA))
          } else {
            grid.text(sprintf("%.1f", cor_geo_bray[i, j]), x, y, gp = gpar(fontsize = 10))
          }
        }, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE)

dev.off()


## PCoA

library(ape)
library(vegan)

# Order the rows and columns
od <- hclust(dist(df.dist.br))$order
m2 <- df.dist.br[od, od]
ordination <- pcoa(m2, correction = "none")
scores <- as.data.frame(ordination$vectors[,1:2]) %>% rownames_to_column(var = "sample.id")
proportion_exp <- ordination$values$Relative_eig[1:2]

#Join metadata to coordinates
atlan.bray <- merge(scores, sample_ste[,c("sample.id","site")], 
                     by = "sample.id", all.x = TRUE)

# Set colors
col = c("#1B9E77","#E6F598","#CA9822", "#351B9E")

pdf("Stereo/PCoA_stereo_atlanticum.pdf",
    width = 8,
    height = 5)

atlan.bray %>% 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(fill = site), size = 3, pch = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dotted") +  # Add vertical line at x = 0
  scale_fill_manual(values = col, name = "Collection site") +
  theme_minimal() +
  scale_size(guide = "none") +  # Remove legend
  labs(x = paste0("Axis 1 (", round(proportion_exp[1] * 100), "%)"),  # Adjusted axis labels with percentages
       y = paste0("Axis 2 (", round(proportion_exp[2] * 100), "%)"),  # Adjusted axis labels with percentages
       title = "PCoA Stereocaulon atlanticum") +  # Improved title
  theme(axis.title = element_text(size = 12),  # Increased size of axis labels
        plot.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Increased size and bold title

dev.off()


