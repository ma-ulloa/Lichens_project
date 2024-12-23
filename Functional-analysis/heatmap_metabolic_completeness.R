library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


#### Heatmap###############################################################

#################################################################################

metabolic_pathways_long <- read.table("kegg-metabolism_modules.txt", header = TRUE, sep = "\t") |>
  mutate(Broader_Group = case_when(
    str_detect(genome_name, "sticta_") ~ "Sticta",
    str_detect(genome_name, "stereocaulon_") ~ "Stereocaulon",
    str_detect(genome_name, "cora_") ~ "Cora",
    TRUE ~ "Other")) |>
    mutate(taxa = sub(".*_", "", genome_name),
    genome_name = gsub("_", " ", genome_name)) |>
    as_tibble()

# Create a mapping for unique Generic_Names
unique_genome_names <- metabolic_pathways_long |>
  distinct(Broader_Group, genome_name) |>
  group_by(Broader_Group) |>
  mutate(Generic_Name = paste(Broader_Group, row_number(), sep = " ")) |>
  ungroup()

# Join the Generic_Names back to the main table
metabolic_pathways_long <- metabolic_pathways_long |>
  left_join(unique_genome_names, by = c("Broader_Group", "genome_name")) |>
  mutate(
    taxa = ifelse(taxa == "pseudomonadota", "Pseudomonadota", taxa)
  ) |> filter(module_category != "Gene set", module_category != "Module set",
              module_category !="Lipid metabolism", module_category != "Nucleotide metabolism",
              module_category != "Glycan metabolism", module_category != "Xenobiotics biodegradation",
              module_category != "Biosynthesis of terpenoids and polyketides", 
              module != "M00882")



module_matrix <- metabolic_pathways_long |>
  select(Generic_Name, module, pathwise_module_completeness) |> # Select relevant columns
  pivot_wider(
    names_from = module,
    values_from = pathwise_module_completeness,
    values_fill = 0 # Fill missing values with 0
  )


##########################

# Set genome_name as rownames
module_matrix <- column_to_rownames(module_matrix, var = "Generic_Name")
module_matrix <- module_matrix |> t()

##########################

anot_df <- metabolic_pathways_long |> distinct(Broader_Group, Generic_Name, taxa) 
col_order <- anot_df$Generic_Name

# Row annotations for MAGs
annotation <- columnAnnotation(
  Taxonomy = anot_df$taxa,
  Lichen = anot_df$Broader_Group,
  col = list(
    Taxonomy = c(
      "Basidiomycota" = "#2A355A", 
      "Burkholderiales" = "#833984",
      "Acetobacteraceae" = "#8E45B0", 
      "Caulobacteraceae" = "#5A2A4D", 
      "Cyanobacteriota" = "#15CB9B", 
      "Beijerinckiaceae" = "#BBC1E2", 
      "Pseudomonadaceae" = "#321A24", 
      "Hyphomicrobiales" = "#B05645", 
      "Burkholderiaceae" = "#845C39",
      "Chlorophyta" = "darkgreen",
      "Acidobacteriota" = "#5A4F2A", 
      "Sphingomonadales" = "#A3058E", 
      "Bacteroidota" = "#ACD392", 
      "Ascomycota" = "#459DB0", 
      "Acidobacteriaceae" = "lavender", 
      "Pseudomonadota" = "gold",
      "Bacteroidota" = "#6e0d25",
      "Armatimonadota" = "#ffffb3",
      "Lichenibacteriaceae" = "#775253",
      "Methylobacteriaceae" = "#ba5a31",
      "Pseudomonas" = "#CF8BA9",
      "Verrucomicrobiota" = "#86e7b8",
      "Sphingomonadaceae" = "cornflowerblue"
    ),
    Lichen = c(
      "Cora" = "#1B9E77", "Stereocaulon" = "#D95F02", "Sticta" = "#7570B3"
    )
  )
)

# Ensure that module_info aligns with the filtered matrix
module_info <- metabolic_pathways_long |> 
  select(module, module_category) |> 
  distinct() 

module_info <- module_info |> filter(module %in% rownames(module_matrix))

# Define colors for each metabolic pathway category
pathway_colors <- c(
  "Carbohydrate metabolism" = "#FFA07A",
  "Energy metabolism" = "#20B2AA",
  "Lipid metabolism" = "#9370DB",
  "Nucleotide metabolism" = "#FFD700",
  "Amino acid metabolism" = "#FF69B4",
  "Metabolism of cofactors and vitamins" = "#8FBC8F",
  "Biosynthesis of terpenoids and polyketides" = "#6495ED",
  "Glycan metabolism" = "#DAA520",
  "Biosynthesis of other secondary metabolites" = "#FF6347",
  "Xenobiotics biodegradation" = "#af1b3f"#,
)

# Create row annotation with the custom colors
row_annotation <- rowAnnotation(
  Category = module_info$module_category,
  col = list(Category = pathway_colors)
)

# Create a grayscale color scale for the heatmap fill
col_fun <- colorRamp2(c(0, 1), c("white", "black"))


# Generate the heatmap
heatmap <- Heatmap(
  module_matrix,
  name = "Metabolic Pathway Completeness",
  col = col_fun,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),  # Set column names font size
  heatmap_legend_param = list(title_gp = gpar(fontsize = 12)),  # Set legend title font size
  column_order = col_order,             # Set the desired fixed column order
  row_title_gp = gpar(fontsize = 0),  # Set row title font size
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_split = anot_df$Broader_Group,
  row_split = module_info$module_category,
  left_annotation = row_annotation,
  top_annotation = annotation,         # Column annotation for modules
  show_row_names = TRUE,
  show_column_names = TRUE
)

# Draw the heatmap
draw(heatmap)


