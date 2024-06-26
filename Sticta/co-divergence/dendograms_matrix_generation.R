---
#title: "Dendograms and Mantel test"
#author: "Alejandra Ulloa"
  
# **1. MANTEL TEST** 

## Mantel tests and mantel correlograms of multivariate dissimilarities vs. phylogenetic distances#


#Load packages
library(correlation)
library(ape)
library(vegan)
library(phangorn)
library(phyloseq)

### Generate ITS phylogeny distance matrix using the cophenetic.phylo function from the ape package

# Load the phylogenetic tree of the Sticta genus
tree.st <- read.tree("Sticta/co-divergence/Reviewed/ht.newick")
plot(tree.st)

ITS_matrix <- cophenetic.phylo(tree.st) #Generate distance matrix from phylogenetic distances

#Load relative abundance tables collapsed at different taxonomic levels
ral.abd.asv <- read.table("Sticta/co-divergence/Reviewed/Abundance_tables/compiled_table_asv.tsv",
                         sep = "\t", header = TRUE, row.names = 1)

ral.abd.gn <- read.table("Sticta/co-divergence/Reviewed/Abundance_tables/compiled_table_genus.tsv",
                         sep = "\t", header = TRUE, row.names = 1)

ral.abd.fm <- read.table("Sticta/co-divergence/Reviewed/Abundance_tables/compiled_table_family.tsv",
                         sep = "\t", header = TRUE, row.names = 1)

ral.abd.cl <- read.table("Sticta/co-divergence/Reviewed/Abundance_tables/compiled_table_class.tsv",
                         sep = "\t", header = TRUE, row.names = 1)

#### Bray-Curtis matrix

#Calculate Bray-Curtis distance for relative abundance table at different tax levels

#Class level
bray_dist_cl <- vegdist(t(ral.abd.cl), method="bray")
bray_dist_mt_cl <- as.matrix(bray_dist_cl) #Generate matrix

#Family level
bray_dist_fm <- vegdist(t(ral.abd.fm), method="bray")
bray_dist_mt_fm <- as.matrix(bray_dist_fm) #Generate matrix

#Genus level
bray_dist_gn <- vegdist(t(ral.abd.gn), method="bray")
bray_dist_mt_gn <- as.matrix(bray_dist_gn) #Generate matrix

#ASV level
bray_dist_asv <- vegdist(t(ral.abd.asv), method="bray")
bray_dist_mt_asv <- as.matrix(bray_dist_asv) #Generate matrix

# Unifrac matrix

#Create phyloseq object with 16tree and ASV relative abundance table of samples

tree_16S = read.tree("Raw_data/tree-16S.nwk")
OTU_16S = otu_table(ral.abd.asv, taxa_are_rows = TRUE)
unifrac.ps <- phyloseq(OTU_16S, tree_16S)
unfrac <- UniFrac(unifrac.ps, FALSE)
weighted_unifrac <- UniFrac(unifrac.ps, TRUE)

# Create Phyloseq objects

####### Calculate similarity of matrices with the mantel() function ########

mantel_bray_cl = vegan::mantel(bray_dist_mt_cl, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
#mantel_bray_cl
#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#vegan::mantel(xdis = bray_dist_mt_cl, ydis = ITS_matrix, method = "pearson",      permutations = 9999, na.rm = TRUE) 

#Mantel statistic r: 0.07635 
#      Significance: 0.31012 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.317 0.475 0.520 0.546 
#Permutation: free
#Number of permutations: 5039

mantel_bray_fm = vegan::mantel(bray_dist_mt_fm, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
#> mantel_bray_fm

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#vegan::mantel(xdis = bray_dist_mt_fm, ydis = ITS_matrix, method = "pearson",      permutations = 9999, na.rm = TRUE) 

#Mantel statistic r: 0.03538 
#      Significance: 0.43492 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.324 0.394 0.472 0.535 
# Permutation: free
# Number of permutations: 5039

mantel_bray_gn = vegan::mantel(bray_dist_mt_gn, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
#> mantel_bray_gn

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#vegan::mantel(xdis = bray_dist_mt_gn, ydis = ITS_matrix, method = "pearson",      permutations = 9999, na.rm = TRUE) 

#Mantel statistic r: -0.008722 
#      Significance: 0.49623 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.319 0.392 0.438 0.481 
#Permutation: free
#Number of permutations: 5039

mantel_bray_asv = vegan::mantel(bray_dist_mt_asv, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
#mantel_bray_asv

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#vegan::mantel(xdis = bray_dist_mt_asv, ydis = ITS_matrix, method = "pearson",      permutations = 9999, na.rm = TRUE) 

#Mantel statistic r: 0.0009607 
#      Significance: 0.46389 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.319 0.372 0.433 0.485 
#Permutation: free
#Number of permutations: 5039

mantel_unifrac_asv = vegan::mantel(unfrac, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel_unifrac_asv = vegan::mantel(weighted_unifrac, ITS_matrix, method = "pearson", permutations = 9999, na.rm = TRUE)


######################## 2. Topological congruence #########################

## Microbiota Dendrograms

###################### Generating Cladogram for the matrices ######################

#Bray-Curtis Class level
px.br.cl <- phangorn::upgma(bray_dist_cl) # UPGMA - Ultrametric 
ape::write.tree(px.br.cl, file = "Sticta/co-divergence/Reviewed/bray-curtis-collapse-class.newick")

#Bray-Curtis Family level
px.br.fm <- phangorn::upgma(bray_dist_fm) # UPGMA - Ultrametric 
ape::write.tree(px.br.fm, file = "Sticta/co-divergence/Reviewed/bray-curtis-collapse-family.newick")

#Bray-Curtis Genus level
px.br.gn <- phangorn::upgma(bray_dist_gn) # UPGMA - Ultrametric  
ape::write.tree(px.br.gn, file = "Sticta/co-divergence/Reviewed/bray-curtis-collapse-genus.newick")

#Bray-Curtis ASV level
px.br.asv <- phangorn::upgma(bray_dist_asv) # UPGMA - Ultrametric  
ape::write.tree(px.br.asv, file = "Sticta/co-divergence/Reviewed/bray-curtis-collapse-asv.newick")

#Unifrac
px.uni.asv <- phangorn::upgma(unfrac) # UPGMA - Ultrametric  
ape::write.tree(px.uni.asv, file = "Sticta/co-divergence/Reviewed/unifrac-collapse-asv.newick")

#Weighted Unifrac
#Unifrac
px.wuni.asv <- phangorn::upgma(weighted_unifrac) # UPGMA - Ultrametric  
ape::write.tree(px.wuni.asv, file = "Sticta/co-divergence/Reviewed/wunifrac-collapse-asv.newick")


