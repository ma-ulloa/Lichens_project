
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
library(metagMisc)
library(phyloseqCompanion)
library(dplyr)

######## Import clean phyloseq object ##########

load("Lichens_clean.RData")


########################## Common and distinctive ASVs #########################
################################################################################

#With Phyloseq object select the ASV relative abundance table
#create and absence/presence table

#venn_ps <- phyloseq_standardize_otu_abundance(ps, method = "pa")

#Taxonomy
tx_ps_new <- taxa.data.table(ps.rarefied)

#Leave only those with a rowmeans > 0
otu_venn <- transform_sample_counts(ps.rarefied, function(x) x / sum(x))
otu_venn <- otu.data.table(otu_venn)
otu_venn <- column_to_rownames(otu_venn, var = "Sample") 
otu_venn <- otu_venn[rowSums(otu_venn[])>0,]
otu_venn <- otu_venn %>% filter_all(any_vars(. > 0.01)) #Select all asvs with a rel abd > 0.01



#Select by genus

otu_venn$Cora <- rowMeans(otu_venn[,c("Sample1", "Sample9","Sample21","Sample23",
                                      "Sample33", "Sample34","Sample35", "Sample47", 
                                      "Sample55","Sample57","Sample72",
                                      "Sample74","Sample75","Sample76","Sample77",
                                      "Sample83","Sample84","Sample116",
                                      "Sample118","Sample119","Sample120","Sample121",
                                      "Sample122","Sample123","Sample124","Sample125",
                                      "Sample163","Sample164","Sample165","Sample166",
                                      "Sample167","Sample168","Sample169","Sample170",
                                      "Sample171","Sample172")], na.rm = TRUE)

otu_venn$Sticta <- rowMeans(otu_venn[,c("Sample4", "Sample8","Sample14","Sample15",
                                        "Sample20", "Sample22","Sample29", "Sample30", 
                                        "Sample31","Sample37","Sample38","Sample39",
                                        "Sample40","Sample41","Sample42","Sample58",
                                        "Sample59","Sample60","Sample61","Sample65",
                                        "Sample66","Sample67","Sample69","Sample70",
                                        "Sample71","Sample101","Sample103","Sample105",
                                        "Sample106","Sample107","Sample108","Sample109",
                                        "Sample110","Sample111","Sample112","Sample113",
                                        "Sample114","Sample115","Sample129","Sample130",
                                        "Sample148","Sample149","Sample150","Sample151",
                                        "Sample153","Sample154","Sample155","Sample156",
                                        "Sample158","Sample159","Sample160","Sample161"
)], na.rm = TRUE)

otu_venn$Stereocaulon <- rowMeans(otu_venn[,c("Sample5", "Sample13","Sample26","Sample27",
                                              "Sample28", "Sample45","Sample49", "Sample51", 
                                              "Sample80","Sample85","Sample86",
                                              "Sample89","Sample90","Sample91","Sample92",
                                              "Sample93","Sample95","Sample96","Sample97",
                                              "Sample99","Sample100","Sample135",
                                              "Sample137","Sample139","Sample140",
                                              "Sample141","Sample142","Sample143","Sample144",
                                              "Sample145" )], na.rm = TRUE)

otu_v <- otu_venn[,c(119:121)]
otu_v2 <- otu_v %>% mutate_if(is.numeric, ~1 * (. > 0))

## Create venn diagram

cora.v <- subset(otu_v2, Cora == 1) 
Co <- row.names(cora.v)

sticta.v <- subset(otu_v2, Sticta == 1)
St <- row.names(sticta.v)

stereo.v <- subset(otu_v2, Stereocaulon == 1)
S <- row.names(stereo.v)

## Create list
venn_lichen <- list(Cora = Co, Sticta = St, Stereocaulon = S)

##Calculate overlap
overlap <- VennDiagram::calculate.overlap(venn_lichen)

#Join areas to taxonomyfile
share.asv <- data.frame(sample = overlap$a5) %>% left_join(tx_ps_new, share.asv, by= "sample")
sticta.unique <- data.frame(sample = overlap$a3) %>% left_join(tx_ps_new, share.asv, by= "sample")
cora.unique <- data.frame(sample = overlap$a1) %>% left_join(tx_ps_new, share.asv, by= "sample")
stereo.unique <- data.frame(sample = overlap$a7) %>% left_join(tx_ps_new, share.asv, by= "sample")


library("ggvenn")
# Default plot

venn.plot <- "Meta-analysis/Venn-plot.pdf"


pdf(file = venn.plot,   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

ggvenn(venn_lichen, 
       fill_color = c("#1B9E77", "#7570B3", "#D95F02", "#CD534CFF"),
       stroke_size = 0.6, set_name_size = 4)

dev.off()
