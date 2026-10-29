##Sangita Abhrajyoti Sir
##Metagenomics data analyis
library(phyloseq)
library(MicrobiotaProcess)
library(ggplot2)
# library(microbiomeutilities)
library(microeco)
library(file2meco)
library(tidyverse)

setwd("/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Wgs/bracken_analysis/")


ps <- import_biom("bracken.biom1", parseFunction = parse_taxonomy_greengenes)

metadata <- read.csv("metadata_bracken.csv", row.names = 1)
sample_data(ps) <- metadata
ps_seed <- subset_samples(ps, Group != "Seed")
ps_seed_control_soil <- subset_samples(ps_seed, Group != "Control_Soil")

# mt <- phyloseq2meco(ps)
# mt$tidy_dataset()
# 
# mt$filter_pollution(taxa = c("Chordata", "Chlorophyta", "Streptophyta"))
# mt$tidy_dataset()


rownames(ps_seed_control_soil@sam_data) <- ps_seed_control_soil@sam_data$Conditions
colnames(ps_seed_control_soil@otu_table) <- ps_seed_control_soil@sam_data$Conditions

# meta <- read.csv("metadata.csv", row.names = 1)
# sample_data(ps_seed_control_soil) <- meta


ps_seed_control_soil <- subset_taxa(ps_seed_control_soil, Phylum != "Chordata")
ps_1 <- subset_taxa(ps_seed_control_soil, Phylum != "Chlorophyta")
ps_2 <- subset_taxa(ps_1, Phylum != "Streptophyta")

rownames(ps_2@sam_data) <- ps_2@sam_data$Conditions
colnames(ps_2@otu_table) <- ps_2@sam_data$Conditions
# 
# meta <- read.csv("metadata.csv", row.names = 1)
# sample_data(ps_seed_control_soil) <- meta


ps_species <- tax_glom(ps_2, "Species", NArm = TRUE)

ps_tax <- as.data.frame(ps_species@tax_table)
ps_tax$Species_1 <- paste(ps_tax$Genus, ps_tax$Species)
ps_tax <- ps_tax[,-7]
colnames(ps_tax)[7] <- "Species"
ps_tax <- as.matrix(ps_tax)
tax_table(ps_species) <- ps_tax



# library(file2meco)
# library(microeco)
# 
# mt <- phyloseq2meco(ps_2)
# 
# t1 <- trans_abund$new(dataset = mt, taxrank = "Species", ntaxa = 25, groupmean = "Conditions")
# g1 <- t1$plot_bar(bar_full = FALSE, legend_text_italic = FALSE)
# g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))





# convert phyloseq object to mpse
mpse4 <- ps_species %>% as.MPSE()

mpse4$Group <- factor(
  mpse4$Group,
  # levels = c("Control_Soil",
  #   "Bulk_Soil",
  #   "Rhizosphere_Soil")
  levels = c("Bulk_Soil",
             "Rhizosphere_Soil")
)


# mpse4 %<>%
#   mp_cal_abundance( # for each samples
#     .abundance = Abundance,
#     force = TRUE
#   ) %>%
#   mp_cal_abundance( # for each groups
#     .abundance=Abundance ,
#     .group=Conditions,
#     force = TRUE
#   )
# mpse4


# mpse_reordered <- mpse4 %>%
#   mutate(
#     Group = factor(
#       Group,
#       levels = c(
#         "0D_Control",
#         "0D_AB211",
#         "0D_AB211_RE",
#         "7D_Control",
#         "7D_AB211",
#         "7D_AB211_RE",
#         "14D_AB211",
#         "14D_AB211_RE",
#         "21D_Control",
#         "21D_AB211",
#         "21D_AB211_RE",
#         "Control_Soil",
#         "Seed"
#       )
#     )
#   )


# visualize the relative abundance of top 10 phyla for each sample.
p1 <- mpse4 %>%
  mp_plot_abundance(
    .abundance=Abundance,
    .group=Group,
    taxa.class = Species,
    topn = 12,
    relative = TRUE,
    force = TRUE,
    order.by.feature = FALSE,
    plot.group = FALSE,
    geom = "flowbar"
  ) +theme(
    # Axis titles and text
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black")) + 
  # guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = "right",
    # legend.justification = "center",
    # legend.box = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  )
p1
ggsave(p1, file = "top_12_species_updated.svg", height = 8, width = 18)

# visualize the relative abundance of top 10 Species of Actinomycetota for each sample.

ps_species_actinomycetota <- subset_taxa(ps_species, Phylum == "Actinomycetota")
# ps_species_actinomycetota <- subset_samples(ps_species_actinomycetota, Conditions != "Control_Soil")
ps_species_actinomycetota <- subset_samples(ps_species_actinomycetota, Conditions != "Seed")

# convert phyloseq object to mpse
mpse4_actinomycetota <- ps_species_actinomycetota %>% as.MPSE()
# visualize the relative abundance of top 10 phyla for each sample.

mpse4_actinomycetota$Group <- factor(
  mpse4_actinomycetota$Group,
  # levels = c("Control_Soil",
  #            "Bulk_Soil",
  #            "Rhizosphere_Soil")
  levels = c("Bulk_Soil",
             "Rhizosphere_Soil")
)


p2 <- mpse4_actinomycetota %>%
  mp_plot_abundance(
    .abundance=Abundance,
    .group=Group,
    taxa.class = Species,
    topn = 10,
    relative = TRUE,
    force = TRUE,
    
  ) +theme(
    # Axis titles and text
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black")) + 
  # guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = "right",
    # legend.justification = "center",
    # legend.box = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  )+ ggtitle("(a.) Top 10 species for phylum Actinomycetota")

p2

ggsave(p2, file = "top_10_Actinomycetota_updated.svg", height = 8, width = 18)

# visualize the relative abundance of top 10 Species of Pseudomonadota for each sample.

ps_species_Pseudomonadota <- subset_taxa(ps_species, Phylum == "Pseudomonadota")
# ps_species_Pseudomonadota <- subset_samples(ps_species_Pseudomonadota, Conditions != "Control_Soil")
ps_species_Pseudomonadota <- subset_samples(ps_species_Pseudomonadota, Conditions != "Seed")

# convert phyloseq object to mpse
mpse4_Pseudomonadota <- ps_species_Pseudomonadota %>% as.MPSE()
# visualize the relative abundance of top 10 phyla for each sample.

mpse4_Pseudomonadota$Group <- factor(
  mpse4_Pseudomonadota$Group,
  # levels = c("Control_Soil",
  #            "Bulk_Soil",
  #            "Rhizosphere_Soil")
  levels = c("Bulk_Soil",
             "Rhizosphere_Soil")
)

p3 <- mpse4_Pseudomonadota %>%
  mp_plot_abundance(
    .abundance=Abundance,
    .group=Group,
    taxa.class = Species,
    topn = 10,
    relative = TRUE,
    force = TRUE,
    
  ) +theme(
    # Axis titles and text
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black")) + 
  # guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = "right",
    # legend.justification = "center",
    # legend.box = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  ) + ggtitle("(a.) Top 10 species for phylum Pseudomonadota")

ggsave(p3, file = "top_10_Pseudomonadota_updated.svg", height = 8, width = 18)

# visualize the relative abundance of top 10 Species of Bacillota for each sample.

ps_species_Bacillota <- subset_taxa(ps_species, Phylum == "Bacillota")
# ps_species_Bacillota <- subset_samples(ps_species_Bacillota, Conditions != "Control_Soil")
ps_species_Bacillota <- subset_samples(ps_species_Bacillota, Conditions != "Seed")

# convert phyloseq object to mpse
mpse4_Bacillota <- ps_species_Bacillota %>% as.MPSE()
# visualize the relative abundance of top 10 phyla for each sample.

mpse4_Bacillota$Group <- factor(
  mpse4_Bacillota$Group,
  levels = c("Bulk_Soil",
             "Rhizosphere_Soil")
)


p4 <- mpse4_Bacillota %>%
  mp_plot_abundance(
    .abundance=Abundance,
    .group=Group,
    taxa.class = Species,
    topn = 10,
    relative = TRUE,
    force = TRUE,
    
  ) +theme(
    # Axis titles and text
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black")) + 
  # guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = "right",
    # legend.justification = "center",
    # legend.box = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  ) + ggtitle("(b.) Top 10 species for phylum Bacillota")

ggsave(p4, file = "top_10_Bacillota_updated.svg", height = 8, width = 18)



