library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

setwd("JOBS/Liquenes/Functional-anlysis/Binning/")

# Metadata

metadata <- read.table("metadata_lichens.txt", header = TRUE, sep = "\t")

#bins <- metadata |> filter(!taxonomical_annotation == "None" & !taxonomical_annotation == "Eukarya")

metabolic_completeness <- read.table("lichens-module_pathwise_completeness-MATRIX.txt", 
                                     header = TRUE, sep = "\t", row.names = 1)

# Remove non present bins 
metabolic_completeness <- metabolic_completeness |> select(-c(stereocaulon_Caulobacterales, stereocaulon_Enterobacterales,
                                                              stereocaulon_Basidiomycota, sticta_Terriglobales, 
                                                              cora_Enterobacterales, cora_Micrococcales, cora_Sphingomonadales))

# Define bin names and corresponding taxa levels
bin_names <- colnames(metabolic_completeness)

# Generate new names by simplifying and numbering
new_bin_names <- c("Cora1", "Cora2", "Cora3", "Cora4", "Cora5", "Cora6", 
                   "Cora7", "Cora8", "Cora9", "Cora10", "Stereo1", 
                   "Stereo2", "Stereo3", "Stereo4","Stereo5", "Stereo6",
                   "Stereo7", "Stereo8", "Sticta1", "Sticta2", "Sticta3", "Sticta4")

# Create a mapping table
mapping_table <- data.frame(
  old_name = bin_names,
  new_name = new_bin_names,
  taxon = gsub(".*_", "", bin_names)  # Extract taxon from original names
)

mapping_table <- left_join(mapping_table, metadata, by = c("old_name" = "name"))

# Rename the columns based on the mapping
colnames(metabolic_completeness) <- mapping_table$new_name[match(colnames(metabolic_completeness), mapping_table$old_name)]

# Enrichement
metabolic_enrichment <- read.table("metabolic_enrichment2.txt", header = TRUE, sep = "\t")

pathways <- read.table("kegg-metabolism_modules.txt", header = TRUE, sep = "\t")
module_info <-  pathways |> select(module, module_category)
module_info <- module_info %>% filter(module %in% metabolic_enrichment$accession) |> unique()

# Keep only enrich modules for heatmap
relevant_ones <- metabolic_enrichment %>% filter(adjusted_q_value < 0.05) %>% select(accession)

#metabolic_landscape <- metabolic_completeness[relevant_ones$accession, bins$name]
metabolic_landscape <- metabolic_completeness[metabolic_enrichment$accession, ]
metabolic_landscape <- as.matrix(metabolic_landscape)


# Ensure metabolic_completeness is a data frame
#metabolic_completeness_df <- metabolic_completeness %>%
#  rownames_to_column("module")  # Convert row names to a column

# Perform the left join to merge metabolic_enrichment with MODULE information
#module_metadata <- metabolic_enrichment %>% filter(accession %in% metabolic_completeness_df$module) |> 
#  select(MODULE, accession)

#safe module metada
#write.table(module_metadata, "module_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

order <- c("Cora5","Cora1", "Cora6", "Cora8", "Cora4", "Cora3", "Cora2", "Cora10","Cora7", "Cora9",
               "Stereo3", "Stereo1","Stereo4","Stereo6","Stereo5","Stereo2","Stereo8","Stereo7","Sticta2",
               "Sticta1", "Sticta4", "Sticta3")

# Set up annotations for the rows (MAGs)
annotation <- columnAnnotation(
  Taxonomy = mapping_table$taxon,
  Lichen = mapping_table$group,
  col = list(
    Taxonomy = c(
      "Basidiomycota" = "#2A355A", 
      "Burkholderiales" = "#833984", 
      "Capsulimonadales" = "#8E45B0", 
      "Caulobacterales" = "#5A2A4D", 
      "Chlorophyta" = "#67B045", 
      "Cyanobacteriota" = "#15CB9B", 
      "Cytophagales" = "#BBC1E2", 
      "Enterobacterales" = "#321A24", 
      "Hyphomicrobiales" = "#B05645", 
      "Micrococcales" = "#845C39", 
      "Rhodospirillales" = "#5A4F2A", 
      "Sphingomonadales" = "#A3058E", 
      "Terriglobales" = "#ACD392", 
      "Ascomycota" = "#459DB0", 
      "Verrucomicrobiota" = "lavender", 
      "Pseudomonadales" = "gold" ),
    Lichen = c(
      "Cora" = "#1B9E77", "Stereocaulon" = "#D95F02", "Sticta" = "#7570B3"
    )
  )
)


# Define colors for each metabolic pathway category
pathway_colors <- c(
  "Carbohydrate metabolism" = "#FFA07A",
  "Energy metabolism" = "#20B2AA",
  "Lipid metabolism" = "#9370DB",
  "Nucleotide metabolism" = "#FFD700",
  "Amino acid metabolism" = "#FF69B4",
  "Metabolism of cofactors and vitamins" = "#8FBC8F",
  "Biosynthesis of terpenoids and polyketides" = "#6495ED",
  "Module set" = "#FF4500",
  "Glycan metabolism" = "#DAA520",
  "Xenobiotics biodegradation" = "#FF6347"
)

# Create row annotation with the custom colors
row_annotation <- rowAnnotation(
  Category = module_info$module_category,
  col = list(Category = pathway_colors)
)

# Create the heatmap
# Create a grayscale color scale for the heatmap fill
col_fun <- colorRamp2(c(0, 1), c("lightgray", "black"))  # Adjust range if needed


# Create the heatmap with grayscale colors and adjusted legend title
heatmap <- Heatmap(
  metabolic_landscape,
  name = "Metabolic Pathway Completeness",
  col = col_fun,
  row_names_gp = gpar(fontsize = 12),
  column_order = order,              # Set the desired fixed column order
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_split = mapping_table$group, # Split columns based on mycobiont groups
  left_annotation = row_annotation,    # Row annotation for taxonomy and mycobiont
  top_annotation = annotation,         # Column annotation for modules
  show_row_names = TRUE,
  show_column_names = TRUE
)

# Draw the heatmap
draw(heatmap)

# Specify the output PDF file name and dimensions for portrait orientation
pdf("heatmap_publication_ready.pdf", width = 8.5, height = 15)  # Adjust width and height as needed for portrait orientation

# Create the heatmap with font size set to 12
heatmap <- Heatmap(
  metabolic_landscape,
  name = "Metabolic Pathway Completeness",
  col = col_fun,
  row_names_gp = gpar(fontsize = 10),  # Set row names font size
  column_names_gp = gpar(fontsize = 10),  # Set column names font size
  heatmap_legend_param = list(title_gp = gpar(fontsize = 12)),  # Set legend title font size
  row_title_gp = gpar(fontsize = 0),  # Set row title font size
  column_title_gp = gpar(fontsize = 12),  # Set column title font size
  column_order = order,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_split = mapping_table$group,      # Cluster rows based on mycobiont groups
  row_split = module_info$module_category,  # Split rows based on module categories
  left_annotation = row_annotation,        # Row annotation for taxonomy and mycobiont
  top_annotation = annotation,  # Column annotation for modules
  show_row_names = TRUE,
  show_column_names = TRUE,
)

# Draw the heatmap
draw(heatmap)

# Close the PDF device
dev.off()


### Compile table

final_table <- left_join(pathways, mapping_table, by = c("db_name" = "old_name")) |>
               select(module, genome_name, new_name, module_category, module_subcategory, module_name,
                      enzymes_unique_to_module) |> filter(module %in% metabolic_enrichment$accession) 

final_table <- final_table |> filter(!is.na(new_name))

write.table(final_table, file = "metabolic_pathways.txt", sep = "\t", col.names = TRUE)


#### New heatmap###############################################################

#################################################################################

metabolic_pathways_long <- read.table("kegg-metabolism_modules_new_dec_18.txt", header = TRUE, sep = "\t") |>
  mutate(Broader_Group = case_when(
    str_detect(genome_name, "sticta_") ~ "Sticta",
    str_detect(genome_name, "stereocaulon_") ~ "Stereocaulon",
    str_detect(genome_name, "cora_") ~ "Cora",
    TRUE ~ "Other"
  )) |>
  mutate(
    taxa = sub(".*_", "", genome_name),
    genome_name = gsub("_", " ", genome_name)
  ) |>
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

# Filter columns based on >50% zeros
#columns_to_keep <- colMeans(module_matrix == 0) <= 0.5
#filtered_module_matrix <- module_matrix[, c(TRUE, columns_to_keep)]

# Convert to matrix
#filtered_module_matrix <- as.matrix(filtered_module_matrix) |> t()

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
  #"Module set" = "#FF4500",
  "Glycan metabolism" = "#DAA520",
  "Biosynthesis of other secondary metabolites" = "#FF6347",
  "Xenobiotics biodegradation" = "#af1b3f"#,
  #"Gene set" = "#a2ad59"
)

# Create row annotation with the custom colors
row_annotation <- rowAnnotation(
  Category = module_info$module_category,
  col = list(Category = pathway_colors)
)

# Create a grayscale color scale for the heatmap fill
col_fun <- colorRamp2(c(0, 1), c("white", "black"))


# Specify the output PDF file name and dimensions for portrait orientation
pdf("heatmap_publication_ready_specific_pathways.pdf", width = 10, height = 25) 
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

# Close the PDF device
dev.off()


###############################

modules <- metabolic_pathways_long |> select(module, module_name, module_category) |> distinct()


### ENRICHMENT

###############################

### modules

modules <- metabolic_pathways_long |> distinct(module, module_category)


###                
metabolic_enrichment <- read.table("metabolic_enrichment_stereo_vs_cora.txt", header = TRUE, sep = "\t") |>
  left_join(modules , by = c("accession" = "module")) |> 
  filter(MODULE != "Heme biosynthesis, animals and fungi, glycine => heme")
metabolic_enrichment <- metabolic_enrichment %>% filter(adjusted_q_value < 0.05) 
metabolic_enrichment <- metabolic_enrichment |> as_tibble()


metabolic_enrichment_long <- metabolic_enrichment |>
  separate_rows(sample_ids, sep = ",") |>
  mutate(Broader_Group = case_when(
      str_detect(sample_ids, "sticta_") ~ "Sticta",
      str_detect(sample_ids, "stereocaulon_") ~ "Stereocaulon",
      str_detect(sample_ids, "cora_") ~ "Cora",
      TRUE ~ "Other")) |> 
  mutate(sample_ids = gsub("^[^_]+_", "", sample_ids))

colors <- c("Cora" = "#1B9E77", "Stereocaulon" = "#D95F02", "Sticta" = "#7570B3")



enrichment_plot <- ggplot(metabolic_enrichment_long, aes(x = Broader_Group, y = MODULE)) +
  geom_point(aes(size = enrichment_score, color = module_category), alpha = 0.7) +
  geom_text(aes(label = sample_ids), hjust = -0.2, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = pathway_colors) +
  theme_bw() +
  labs(
    title = "Metabolic Pathway Enrichment and Bin Contributions",
    x = "Group",
    y = "Metabolic Pathway",
    size = "Enrichment Score",
    color = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("enrichment_plot_stereo_vs_cora.pdf", enrichment_plot, width = 17, height = 12)

#### Enrichment withing each individual lichen genra

metabolic_enrichment_cora <- read.table("metabolic_enrichment_cora.txt", header = TRUE, sep = "\t") |>
  left_join(modules , by = c("accession" = "module")) |> 
  filter(MODULE != "Heme biosynthesis, animals and fungi, glycine => heme")
metabolic_enrichment_stereo <- read.table("metabolic_enrichment_stereo_alg.txt", header = TRUE, sep = "\t") |>
  left_join(modules , by = c("accession" = "module")) |> 
  filter(MODULE != "Heme biosynthesis, animals and fungi, glycine => heme")
metabolic_enrichment_sticta <- read.table("metabolic_enrichment_sticta.txt", header = TRUE, sep = "\t") |>
  left_join(modules , by = c("accession" = "module")) |> 
  filter(MODULE != "Heme biosynthesis, animals and fungi, glycine => heme")

color_palete <- c("Cyanobacteria" = "#63264a", "Bacteria" = "#fed766", "Fungi" = "#06ba63")

enrichment_plot_within_group <- function(enrichment_table, color_palete, adjusted = TRUE) {
  # Filter based on adjusted q-value if requested
  if (adjusted) {
    enrichment_table <- enrichment_table %>% filter(adjusted_q_value < 0.05)
  }
  
  # Ensure input is a tibble
  enrichment_table <- enrichment_table |> as_tibble()
  
  # Transform data into long format and classify groups
  enrichment_table_long <- enrichment_table |>
    separate_rows(sample_ids, sep = ",") |>
    mutate(Broader_Group = case_when(
      str_detect(sample_ids, "Basidiomycota") | str_detect(sample_ids, "Ascomycota") ~ "Fungi",
      str_detect(sample_ids, "Cyanobacteriota") ~ "Cyanobacteria",
      str_detect(sample_ids, "Chlorophyta") ~ "Algae",
      TRUE ~ "Bacteria"
    )) |> 
    mutate(sample_ids = gsub("^[^_]+_", "", sample_ids)) 
  
  # Create the ggplot
  ggplot(enrichment_table_long, aes(x = Broader_Group, y = MODULE)) +
    geom_point(aes(size = enrichment_score, color = module_category), alpha = 0.7) +
    geom_text(aes(label = sample_ids), hjust = -0.2, size = 3, check_overlap = TRUE) +
    scale_color_manual(values = color_palete) +
    theme_bw() +
    labs(
      title = "Metabolic Pathway Enrichment",
      x = "Group",
      y = "Metabolic Pathway",
      size = "Enrichment Score",
      color = "Metabolic Category"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


enr_cora <- enrichment_plot_within_group(metabolic_enrichment_cora, pathway_colors)
enr_stereo <- enrichment_plot_within_group(metabolic_enrichment_stereo, pathway_colors)
enr_sticta <- enrichment_plot_within_group(metabolic_enrichment_sticta, pathway_colors, adjusted = FALSE)

ggsave("enrichment_within_lichen_cora.pdf", enr_cora, width = 14, height = 11)
ggsave("enrichment_within_lichen_stereo.pdf", enr_stereo, width = 13, height = 11)
ggsave("enrichment_within_lichen_sticta.pdf", enr_sticta, width = 14, height = 11)

#################################################################################

########################### Binning stats ######################################

library(jsonlite)

# Read the data
library(tidyverse)
library(jsonlite)

data <- read_table("bin_stats.analyze.tsv")
#data <- data |> filter(! c("X2","X4","X5","X7","X8","X10", 
#                           "X11","X12", "X14","X15","X17", "X18"))
parsed_data <- data %>%
  mutate(
    Details = gsub("([a-zA-Z0-9_# ]+):", '"\\1":', Details), # Add quotes around keys
    Details = map(Details, ~ fromJSON(.x)) # Parse the JSON strings
  ) %>%
  unnest_wider(Details) # Expand the parsed list into columns

# Function to convert the Details column to a proper dataframe
parse_details <- function(details) {
  # Remove curly braces and split key-value pairs
  details <- gsub("[{}]", "", details)  # Remove the curly braces
  key_value_pairs <- str_split(details, ",\\s*")[[1]]  # Split by commas
  
  # Create a named list of key-value pairs
  parsed_details <- map(key_value_pairs, ~ str_split(.x, ":\\s*")[[1]])
  
  # Convert the list to a data frame
  parsed_df <- as.data.frame(t(do.call(rbind, parsed_details)), stringsAsFactors = FALSE)
  
  # Set column names (the keys from the details)
  colnames(parsed_df) <- parsed_df[1, ]
  parsed_df <- parsed_df[-1, ]  # Remove the first row which contains the keys
  
  # Convert the data to the correct type (numeric where applicable)
  parsed_df <- parsed_df %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)), .names = "{.col}_num"))
  
  return(parsed_df)
}

# Apply the function and merge with the original sample column
data_parsed <- data %>%
  mutate(Details = map(Details, parse_details)) %>%
  unnest(Details)

# View the result
print(data_parsed)

