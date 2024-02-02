
#title: "Liquenes_análisis_16S"
#author: "Alejandra Ulloa"
#date: ' December 2023'

#This script creates the ASV heatmap for the most abundant families and those
#present in the biplot.


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


######################## Heatmap of most abundant ASVs within 8 families #################################

#Calculate the relative abundance
abundant.asv <-  ps.rarefied %>%
  tax_mutate(Species = NULL) %>%
  tax_select(c("f__Acetobacteraceae",
               "f__Acidobacteriaceae_(Subgroup_1)",
               "f__Beijerinckiaceae",
               "f__Caulobacteraceae" ,
               "f__Sphingobacteriaceae",
               "f__Sphingomonadaceae",
               "f__Nostocaceae",
               "f__Deinococcaceae")) %>%
  tax_top(n = 30, by = mean, rank = "unique", use_counts = FALSE)

ps.rarefied.abd <-  ps.rarefied %>%
  tax_mutate(Species = NULL) %>%
  tax_select(abundant.asv) %>%
  tax_fix() %>%
  tax_fix(unknowns = c("f__uncultured")) %>%
  tax_transform("compositional")


##Heatmap with Complex heatmaps
library(ComplexHeatmap)

#Abundance table
abd.tab <- otu.data.table(ps.rarefied.abd) 

# Create a common identifier column
abd.tab$CommonID <- 1:nrow(abd.tab)  
#Add column with new names
abd.tab$ASV_ID <- paste("asv", 1:nrow(abd.tab), sep = "")

#Tax_table add new names to they match with the otu table
tax.t <- taxa.data.table(ps.rarefied.abd)
tax.t$CommonID <- 1:nrow(tax.t)
tax.t$ASV_ID <- paste("asv", tax.t$CommonID, sep = "")

#Edit abundance table
abd.tab <- column_to_rownames(abd.tab, var = "ASV_ID")
abd.tab <- abd.tab[,c(-1,-120)]  #Remove "Sample" and "common identifier" columns
abd.matrix <- as.matrix(abd.tab)

#Heatmap annotations
genus <- sample.data.frame(ps.rarefied.abd)

# Abundance colors
palette_blues <- colorRampPalette(colors = c("white", "#00838EFF"))(9)

#Anotate samples according to the genus they belong to and selection of colors to use
anot = HeatmapAnnotation(Genus = genus$host_genus, col = list(Genus = c("Cora" = "#1B9E77", 
                                                                        "Stereocaulon" = "#D95F02", 
                                                                        "Sticta"= "#7570B3"), which="col"))

#Anotate samples according to the genus they belong to and selection of colors to use
anot2 = HeatmapAnnotation(Species = genus$host_species_ITS, 
                          col = list(Species = c("Cora_paraciferrii" ="#1B9E77",
                                                 "Cora_sp" = "#B2DFDAFF",
                                                 "Cora_celestinoa" = "#7FCBC4FF",
                                                 "Cora_davidia" = "#1B849E",
                                                 "Cora_putumayensis" = "#00685BFF",
                                                 "Stereocaulon_atlanticum" = "#FF9800FF",
                                                 "Stereocaulon_sp" = "#E55100FF",
                                                 "Sticta_arachnofuliginosa" = "#4A138CFF",
                                                 "Sticta_gyalocarpa" = "#8D24AAFF",
                                                 "Sticta_macrofuliginosa" = "#AB46BBFF",
                                                 "Sticta_minutula" = "#B967C7FF",
                                                 "Sticta_rhizinata" = "#CD92D8FF",
                                                 "Sticta_sp"= "#E0BEE6FF",
                                                 "Sticta_sylvatica" = "#F2E5F4FF"), which="col"))

#Anota ASVs according to the family they belong to and selection of colors to use
left.anot = rowAnnotation(Family = tax.t$Family, 
                          col = list(Family = c("f__Acetobacteraceae" = "#351B9E",
                                                "f__Acidobacteriaceae_(Subgroup_1)" = "#B185F0",
                                                "f__Beijerinckiaceae" = "#858FF0",
                                                "f__Deinococcaceae"  = "#3C5488B2",
                                                "f__Caulobacteraceae" = "#F0E685FF",
                                                "f__Sphingobacteriaceae" = "#CA9822",
                                                "f__Sphingomonadaceae" = "#9E771B",
                                                "f__Nostocaceae"= "#849E1B"
                          )))

### Column order
library("seriation")

genus.ord <- genus[order(genus$host_species_ITS),]
genus.col.ord <- genus.ord$sample.id

#row order
asv.order <- tax.t[order(tax.t$Family, decreasing = TRUE),]
asv.row.ord <- asv.order$ASV_ID

#Heatmap

outfile_HM <- "Meta-analysis/Heatmap_ASV.pdf"

pdf(file = outfile_HM ,   # The directory you want to save the file in
    width = 20 ,          # The width of the plot in inches
    height = 4.5,
    colormodel = "cmyk")  # Color model (cmyk is required for most publications

Heatmap(abd.matrix, 
        name = "Relative abundance",
        row_names_side = "left", 
        column_names_side = "bottom", 
        row_dend_side = "left",
        #top_annotation = anot,
        bottom_annotation = anot2,
        left_annotation = left.anot,
        rect_gp = gpar(col = "grey"), 
        col = palette_blues,
        show_row_dend = FALSE,
        row_names_gp = gpar(cex=0.9, fontface = "italic"), 
        column_names_gp = gpar(cex=1), 
        column_dend_height = unit(0.5, "cm"), 
        column_title = "Samples",
        column_names_rot = 90, 
        #clustering_distance_columns = bray_dist, #Insert fuction for clustering
        #row_title = "ASV", 
        #column_order = col.order.genus,
        column_order = genus.col.ord,
        row_order = asv.row.ord)
#clustering_distance_columns = bray_dist,
#clustering_method_rows = "ward.D2"

dev.off()




