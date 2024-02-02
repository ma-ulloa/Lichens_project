
# Title: "Tanglegrams visualization"
# author: Maria Alejandra Ulloa

#Load packages
library(correlation)
library(ape)
library(vegan)
library(phyloseq)
library(BiocManager)
library(dendextend)
library(stats)
library(phangorn)
library(phytools)
library(colorspace)


# 1. MANTEL TEST

# Generate ITS phylogeny distance matrix using 
# cophenetic.phylo function from the ape package

# Load the phylogenetic tree of the Sticta genus
tree.st <- read.tree("Sticta/co-divergence/Sticta-phylogeny/ht.newick")
plot(tree.st)

# Generate distance matrix from phylogenetic distances
ITS_matrix <- cophenetic.phylo(tree.st) 

# Bray-Curtis matrix

ral.abd.asv <- read.table("Sticta/co-divergence/Relative-abundance-tables/Ral-sticta-asv-by-specie.txt")
bray_dist_asv <- vegdist(t(ral.abd.asv), method="bray")
bray_dist_mt_asv <- as.matrix(bray_dist_asv) 

# Mantel
mantel_bray_asv = vegan::mantel(bray_dist_mt_asv, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel_bray_asv

# 2. Topological congruence

# Microbiota Dendrograms

#Bray-Curtis Class level
px.br.cl <- read.tree("Sticta/co-divergence/Dendograms/bray-curtis-collapse-class.newick")

#Bray-Curtis Genus level
px.br.gn <- read.tree("Sticta/co-divergence/Dendograms/bray-curtis-collapse-genus.newick") 

#Bray-Curtis asv level
px.br.asv <- read.tree("Sticta/co-divergence/Dendograms/bray-curtis-collapse-asv.newick") 
plot(px.br.asv)

#Unifrac
px.uni <- read.tree("Sticta/co-divergence/Dendograms/unifrac.newick")
px.wuni <- read.tree("Sticta/co-divergence/Dendograms/wunifrac.newick")

# Make a list with all dendrograms they can be use with the function tanglegram()
microbiota.dend <- c(px.br.asv, px.br.gn, px.uni, px.wuni)

##### Visualization comparison

# ITS Tree pre-processing

# Make tree ultrametric
tpxo <- chronos(tree.st)

# Test wether the tree is ultra metric, binary and rooted
is.ultrametric(tpxo)
is.binary(tpxo)
is.rooted(tpxo)

#Make tree a hierarchical objecct and dendogram object
treeclust <- as.hclust.phylo(tpxo) %>% as.dendrogram()

# save tree as dendongram object
pxo1 <- as.dendrogram(treeclust)
plot(pxo1)


### Loop tanglegrams visualization

pdf("Sticta/co-divergence/tanglegram_Unifrac.pdf", 
    width = 7 , 
    height = 5)

# This function takes as input the microbiota dendrogram and generates a 
# tanglegram comparing the topologies with the mycobiont phylogeny.

tanglegrams = function(tree.2) {
  
  #Make tree a hierarchical objecct and dendrogram
  dend.tree <- tree.2 %>% as.hclust.phylo() %>% as.dendrogram()
  #Intersect trees
  dend_intersect <- intersect_trees(pxo1, dend.tree)
  ## Save pruned tree_1 
  PXO1 <- dend_intersect[[1]]
  ## save pruned tree_2
  PXO.tree <- dend_intersect[[2]]
  ## Match order by labels
  match_order_by_labels(PXO1,PXO.tree)
  ## generates the tanglegram
  some_colors <- rainbow_hcl(7)
  d1_col <- some_colors[order.dendrogram(PXO1)]
  d2_col <- some_colors[order.dendrogram(PXO.tree)]
  
  dl <- dendlist(
    PXO1 %>%
      set("branches_lty", 1) %>%
      set("leaves_pch", 20) %>%
      set("labels_cex", 1) %>%
      set("branches_lwd", 2),
    
    PXO.tree %>% 
      set("branches_lty", 1) %>%
      set("leaves_pch", 20) %>%
      set("labels_cex", 1) %>%
      set("branches_lwd", 2)
  ) 
  
  dendextend::untangle(PXO1, PXO.tree, method = "random") %>% #untangle
    
  #Generate tanglegram
  tanglegram( common_subtrees_color_lines = TRUE, 
              highlight_distinct_edges  = FALSE, 
              highlight_branches_lwd=FALSE, 
              margin_inner=10,
              lwd=2 )
  
}

tanglegrams(px.uni)

#res <- lapply(microbiota.dend, tanglegrams)

dev.off()




