
#title: "Functional-analysis"
#author: "Maria Alejandra Ulloa"
#---

#Load packages

library(data.table)
library(tidyr)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(stringr)


## COG CATEGORIES ANALYSIS

# Load data

cog.data <- read.csv(file = "Functional-analysis/COG_metabolism_frequency.csv", sep = ";")
cog.data <- cog.data[,-23] # Drop unknown cog category column
cog.dic <- readxl::read_xlsx("Functional-analysis/Symbol_cog.xlsx") 
cog.dic <- cog.dic[,c(-3,-4)]

cog.data.long <- cog.data %>% pivot_longer(cols = -Genero,
                                           names_to = "COG Category",
                                           values_to = "No. of genes")
cog.data.long <- left_join(cog.dic, cog.data.long, by = "COG Category")

# Replace dots with spaces in COG Category names
cog.data.long$`COG Category` <- str_replace_all(cog.data.long$`COG Category`, "\\.", " ")

# Calculate total number of genes for each genus and sort by descending order
genus_order <- cog.data.long %>% 
  group_by(Genero) %>% 
  summarize(total_genes = sum(`No. of genes`)) %>% 
  arrange(desc(total_genes)) %>% 
  pull(Genero)

# Convert Genus to a factor with ordered levels
cog.data.long$Genero <- factor(cog.data.long$Genero, levels = genus_order)

# Calculate total number of genes for each COG category and sort by descending order
cog_category_order <- cog.data.long %>%
  group_by(`COG Category`) %>%
  summarize(total_genes = sum(`No. of genes`)) %>%
  arrange(desc(total_genes)) %>%
  pull(`COG Category`)

# Convert COG Category to a factor with ordered levels
cog.data.long$`COG Category` <- factor(cog.data.long$`COG Category`, levels = cog_category_order)

# Plot colors
col.cog <- c("#383B73", "#7176D9", "#627838", "#C2D971",
             "#8B6D2F", "#D9A84E", "#823B39","#1B9E77", "#7B4173",
             "#CE6DBE", "#D95F02","#7570B3", "#E7298A",
             "#66A61E", "#E6AB02", "#A6761D", "#666666",
             "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
             "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
             "#BC80BD", "#CCEBC5", "#FFED6F",
             "#BF5B17"
)


outfile_cog <- "Functional-analysis/plot_cog.pdf"

pdf(file = outfile_cog,   # The directory you want to save the file in
    width = 16,           # The width of the plot in inches
    height = 5,           # Adjusted the height for better readability
    colormodel = "cmyk")  # Color model (cmyk is required for most publications)

cog.data.long %>% 
  ggplot(aes(x = Genero,
             y = `No. of genes`,
             fill = `COG Category`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col.cog) + 
  coord_flip() +                          # Flip the coordinates for a horizontal bar plot
  ggtitle("COG categories associated to genes") + 
  xlab("Genus") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()

dev.off()


