library(tidyverse)
library(patchwork)

setwd("/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Wgs/Functional_analysis/")
process_heatmap <- function(file, top_n = 10, color_high = "red", title_text = "", show_x_labels = TRUE) {
  
  # Load data
  df <- read.csv(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  df <- df[,-17]
  # View(df)
  # Remove Unclassified
  df <- df %>% filter(Species != "Unclassified")
  
  # Reshape to long format
  df_long <- df %>%
    pivot_longer(
      cols = -c(GO_ID, Function, Species),
      names_to = "Sample",
      values_to = "Abundance"
    )
  
  # Select top N species by total abundance
  top_species <- df_long %>%
    group_by(Species) %>%
    summarise(total = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(total)) %>%
    slice(1:top_n) %>%
    pull(Species)
  
  df_top <- df_long %>%
    filter(Species %in% top_species)
  
  # Order samples
  df_top$Sample <- factor(
    df_top$Sample,
    levels = c(
      # "Control_Soil",
      "0D_Control","0D_AB211","0D_AB211_RE",
      "7D_Control","7D_AB211","7D_AB211_RE",
      "14D_Control","14D_AB211","14D_AB211_RE",
      "21D_Control","21D_AB211","21D_AB211_RE"
    )
  )
  
  # Create heatmap
  p <- ggplot(df_top, aes(x = Sample, y = Species, fill = Abundance)) +
    geom_tile(color = "grey40", size = 0.25) +
    scale_fill_gradient(low = "white", high = color_high, name = "CPM", limits = c(0, 5),
                        oob = scales::squish) +
    theme_light(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
      axis.text.x.bottom = if (!show_x_labels) element_blank() else element_text(angle = 45, hjust = 1, face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 10),
      panel.grid = element_blank()
    ) +
    labs(
      x = "",
      y = ""
    ) +
    ggtitle(title_text)
  
  return(p)
}


p1 <- process_heatmap("Nitrogen fixation.csv", top_n = 10, color_high = "blue",
                      title_text = "(a.) Nitrogen fixation (GO:0009399)", 
                      show_x_labels = F)
p1
p2 <- process_heatmap("Carbon fixation.csv", top_n = 10, color_high = "blue",
                      title_text = "(b.) Carbon fixation (GO:0015977)",
                      show_x_labels = F)
p2
p3 <- process_heatmap("BioFilm associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(c.) BioFilm associated functions (GO:0043708,GO:0043709,GO:0044010,GO:0044011,GO:1900189,
                      GO:1900191,GO:1900192,GO:1900232,GO:1900233)",
                      show_x_labels = FALSE)
p3
p4 <- process_heatmap("Chemotaxis associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(d.) Chemotaxis associated functions (GO:0006935, GO:0050920, GO:0050921, GO:0050918)",
                      show_x_labels = F)
p5 <- process_heatmap("Motility associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(e.) Motility associated functions (GO:0043107,GO:0048870,GO:0060296,GO:0071945,GO:0071973,GO:0071977,
                      GO:0071978,GO:0097588,GO:1902093,GO:1902201,GO:2000145,GO:2000147)",
                      show_x_labels = FALSE)
p6 <- process_heatmap("Quorum sensing.csv", top_n = 10, color_high = "blue",
                      title_text = "(f.) Quorum sensing (GO:0009372)",
                      show_x_labels = T)

p6
p7 <- process_heatmap("Tryptophan associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(g.) Tryptophan associated functions (GO:0000162,GO:0004510,GO:0004830,GO:0004833,GO:0004834,GO:0005300,
                      GO:0006436,GO:0006568,GO:0006569,GO:0009034,GO:0015196,GO:0015827,GO:0019441,GO:0019442,
                      GO:0034354,GO:0050361,GO:0052684)",
                      show_x_labels = FALSE)
p8 <- process_heatmap("Siderphore associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(h.) Siderphore associated functions (GO:0015343,GO:0015344,GO:0015891,GO:0019290,GO:1900705)",
                      show_x_labels = FALSE)
p9 <- process_heatmap("ACC Deaminase associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(i.) ACC Deaminase associated functions (GO:0008660)",
                      show_x_labels = FALSE)
p10 <- process_heatmap("Metal ion associated functions.csv", top_n = 10, color_high = "blue",
                      title_text = "(j.) Metal ion associated functions (GO_0000041,GO_0010038,GO_0030001,GO_0046872,GO_0046873,GO_0046914,
                      GO_0051540,GO_0071248,GO_0090423)",
                      show_x_labels = TRUE)

# Combine heatmaps
final_plot <- p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8/ p9/ p10
# final_plot
ggsave(final_plot, filename = "combined_heatmap_updated.svg", height = 22, width = 11, dpi = 600)

final_plot1 <- p1 / p2 
# Save
ggsave(final_plot1, filename = "combined_heatmap_fixation_updated.svg", height = 8, width = 10, dpi = 600)

final_plot2 <- p3 / p4 / p5 / p6
ggsave(final_plot2, filename = "combined_heatmap_motility_updated.svg", height = 13, width = 13, dpi = 600)

final_plot3 <- p7 / p8/ p9/ p10
ggsave(final_plot3, filename = "combined_heatmap_others_updated.svg", height = 13, width = 13, dpi = 600)





#Bacillus aryabhattai
# library(tidyverse)

df <- read.csv("bacillus_aryabhattai.csv", header = TRUE, check.names = FALSE, 
               stringsAsFactors = FALSE, sep = ",")
df <- df [,-18]
# Reshape to long format
df_long <- df %>%
  pivot_longer(
    cols = -c(GO_ID, Function, Label, Species),
    names_to = "Sample",
    values_to = "Abundance"
  )

df_long$Sample <- factor(
  df_long$Sample,
  levels = c(
    "0D_Control","0D_AB211","0D_AB211_RE",
    "7D_Control","7D_AB211","7D_AB211_RE",
    "14D_Control","14D_AB211","14D_AB211_RE",
    "21D_Control","21D_AB211","21D_AB211_RE"
    # "Control_Soil"
  )
)

p <- ggplot(df_long, aes(x = Sample, y = Label, fill = Abundance)) +
  geom_tile(color = "grey40", size = 0.25) +
  scale_fill_gradient(low = "white", high = "red4", name = "CPM", limits = c(0, 5),
                      oob = scales::squish) +
  theme_light(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    # axis.text.x.bottom = if (!show_x_labels) element_blank() else element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = ""
  ) +
  ggtitle("(a.) Gene ontology functions of Bacillus aryabhattai")

p

ggsave(p, filename = "Gene ontology functions of Bacillus aryabhattai 1.svg", height = 8, width = 12, dpi = 600)
























