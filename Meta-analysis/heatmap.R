
#title: "Liquenes_análisis_16S"
#author: "Alejandra Ulloa"


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
anot2 = HeatmapAnnotation(Species = genus$species.checked, 
                          col = list(Species = c("Cora accipiter" ="#00685BFF",
                                                 "Cora aturuoca" = "#1B849E",
                                                 "Cora celestinoa" = "#1B9E77", 
                                                 "Cora davidia" = "#7FCBC4FF",
                                                 "Cora dewisanti" = "#B2DFDAFF",
                                                 "Cora paraciferrii" = "#32533d", 
                                                 "Cora spec Moncada 021" = "#307351",
                                                 "Cora spec Moncada 024" = "#83b692",
                                                 "Cora spec Moncada 095-096" = "#8cc7a1",
                                                 "Cora subdavidcrinita" = "#ccddb7", 
                                                 "Stereocaulon alpinum" = "#E55100FF",
                                                 "Stereocaulon ramulosum" = "#FF9800FF",
                                                 "Stereocaulon novogranatense" = "#fff275",
                                                 "Sticta aff cellulosa" = "#4b2142",
                                                 "Sticta aff gyalocarpa" = "#74226c",
                                                 "Sticta arachnofuliginosa" = "#8D24AAFF",
                                                 "Sticta brevior" = "#AB46BBFF",
                                                 "Sticta rhizinata" = "#4A138CFF",
                                                 "Sticta andina"= "#4F359B",
                                                 "Sticta sylvatica" = "#CD92D8FF",
                                                 "Sticta parahumboldtii" = "#E0BEE6FF",
                                                 "Sticta pseudohumboldtii" ="#F2E5F4FF"), which="col"))

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

genus.ord <- genus[order(genus$species.checked),]
genus.col.ord <- genus.ord$sample.id

#row order
asv.order <- tax.t[order(tax.t$Family, decreasing = TRUE),]
asv.row.ord <- asv.order$ASV_ID

#Heatmap

outfile_HM <- "Meta-analysis/Heatmap_ASV_updated.pdf"

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
        #row_title = "ASV", 
        column_order = genus.col.ord,
        row_order = asv.row.ord)

dev.off()




