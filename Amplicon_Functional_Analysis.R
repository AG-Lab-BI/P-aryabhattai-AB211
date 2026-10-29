library(tidyverse)
library(ggpicrust2)
library(tibble)

# 1. Load data

abundance_file <- "picrust2_input_without_seed.tsv"

metacyc_abundance <- read_delim(
  abundance_file,
  delim = "\t",
  col_names = TRUE,
  trim_ws = TRUE
)
colnames(metacyc_abundance) [1] <- "pathway"

metadata <- read_delim(
  "metadata.csv",
  delim = ",",
  trim_ws = TRUE
)

# 2. Differential abundance analysis

metacyc_daa_results_df <- pathway_daa(
  abundance = metacyc_abundance %>%
    column_to_rownames("pathway"),
  metadata = metadata,
  group = "Group",
  daa_method = "ALDEx2",
  p.adjust = "fdr"
)

annotated_metacyc_daa_results_df <- pathway_annotation(
  pathway = "KO",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = TRUE
)

feature_with_p_0.05 <- annotated_metacyc_daa_results_df %>%
  filter(p_adjust < 0.05)

# 3. Prepare pathway labels

pathway_labels <- feature_with_p_0.05 %>%
  select(feature, pathway_description) %>%
  distinct() %>%
  rename(pathway = feature)

# pathway_labels <- annotated_metacyc_daa_results_df %>%
#   filter(!is.na(p_adjust)) %>%
#   arrange(p_adjust) %>%        # smallest p_adjust first
#   slice_head(n = 50) %>%       # top 50
#   select(feature, description) %>%
#   distinct() %>%
#   rename(pathway = feature)

# ===============================
# 4. Build heatmap abundance matrix
# ===============================
## For metacyc
# df <- pathway_labels %>%
#   inner_join(metacyc_abundance, by = "pathway") %>%
#   column_to_rownames("pathway_description") %>%
#   select(-pathway)

#For Kegg
df <- pathway_labels %>%
  inner_join(metacyc_abundance, by = "pathway") %>%
  filter(!is.na(pathway_description)) %>%
  group_by(pathway_description) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("pathway_description")

write.csv(df, file = "significant_ko.csv")


# 5. Draw heatmap
library(stringr)
library(dplyr)


df <- read.csv("/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Picrust2/selected_pathways_sum_without_seed_control_soil.csv", header = TRUE, sep = ",", row.names = 1)
df <- df[,-40:-42]
df <- as.matrix(df)
mode(df) <- "numeric"
metadata <- read_delim(
  "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Picrust2/metadata.csv",
  delim = ",",
  trim_ws = TRUE
)



# Filter out rows containing "Control_soil" and seed
metadata_clean <- metadata %>% 
  filter(!str_detect(Group, "Control_Soil|Seed"))


# Fix group naming consistency
metadata_clean <- metadata_clean %>%
  mutate(
    Group = gsub("^00D_", "0D_", Group),
    Group = gsub("^07D_", "7D_", Group)
  )

# Desired order
group_order <- c(
  # "Control_Soil",
  
  "0D_Control",
  "0D_AB211",
  "0D_AB211RE",
  
  "7D_Control",
  "7D_AB211",
  "7D_AB211RE",
  
  "14D_Control",
  "14D_AB211",
  "14D_AB211RE",
  
  "21D_Control",
  "21D_AB211",
  "21D_AB211RE"
)

# Set factor levels
metadata_clean$Group <- factor(metadata_clean$Group, levels = group_order)

# Reorder metadata by factor levels
metadata_clean <- metadata_clean %>%
  arrange(Group)

p <- pathway_heatmap(
  abundance = df,
  metadata = metadata_clean,
  group = "Group",
  cluster_rows = F,
  cluster_cols = F,
  clustering_method = "complete",
  clustering_distance = "correlation",
  low_color = "#2166ac",
  mid_color = "#f7f7f7",
  high_color = "#b2182b",
  colorbar_title = "Standardized Abundance",
  font_size = 14)

p

ggsave(p, file = "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Picrust2/picrust2_heatmap_updated_without_seed_control_soil.svg", height = 6, width = 24)



