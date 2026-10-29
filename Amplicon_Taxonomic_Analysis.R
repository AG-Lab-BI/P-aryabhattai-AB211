library(phyloseq)
# library(microbiomeutilities)
library(file2meco)
library(microeco)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(MiscMetabar)
library(ggplot2)
library(tidyverse)
library(treeio)
library(trac)
library(ggnewscale)
library(RColorBrewer)
library(scales)
library(colorspace)


# setwd("D:/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_SIr/Final_tree_phylum")

# ?qiime2meco
otuqzafile <-"/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Qiime2_Outputs/tableNoFilt.qza"
taxaqzafile <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Qiime2_Outputs/taxonomyNoFilt.qza"
mapfile <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Qiime2_Outputs/metadata.tsv"
qiime2meco(otuqzafile)
qiime2meco(otuqzafile, sample_table = mapfile, taxonomy_table = taxaqzafile)
tree_data <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Qiime2_Outputs/treeNoFilt.qza"
rep_data <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Qiime2_Outputs/repseqsNoFilt.qza"


test1 <- qiime2meco(otuqzafile, sample_table = mapfile, taxonomy_table = taxaqzafile, phylo_tree = tree_data, rep_fasta = rep_data, auto_tidy = TRUE)
test1
# ps_tree@tax_table@.Data = substring(ps_tree@tax_table@.Data, 4)

suppressWarnings(suppressMessages(ps <- meco2phyloseq(test1)))
ps@tax_table@.Data = substring(ps@tax_table@.Data, 4)
ps_seed <- subset_samples(ps, Group != "Seed")
ps_seed_cs <- subset_samples(ps_seed, Group != "Control_Soil")

physeq.all <- prune_taxa(taxa_sums(ps_seed_cs)>0, ps_seed_cs) # remove ASVs that are not present anymore
physeq.all_1 <- subset_taxa(physeq.all, Kingdom != "Eukaryota")
physeq.all_2 <- subset_taxa(physeq.all_1, Family != "Mitochondria")
physeq.all_3 <- subset_taxa(physeq.all_2, Order != "Chloroplast")


cat("Nb of Sputum samples:", nsamples(ps))


#_________________________________
## 3.2. Get a treedata object ####
phyloseq_to_treedata <- function(physeq){
  
  ## ***********************
  ## 1 - GET TAXONOMIC TABLE
  tax <- physeq@tax_table@.Data[,1:6] # don't include Species
  # Add a "Life" column (rank0)
  tax <- cbind("Life", tax)
  colnames(tax)[1] <- "Life"
  # tax[1:5,1:7]
  
  
  ## ***********************
  ## 2 - GET TAXONOMIC TABLE WITH CUMULATIVE LABELS
  full_tax <- tax
  for (i in 2:7) {
    full_tax[, i] <- paste(full_tax[, i-1], full_tax[, i], sep = "::")
  }
  full_tax <- as.data.frame(full_tax, stringsAsFactors = TRUE)
  # full_tax[1:5,1:7]
  
  
  ## ***********************
  ## 3 - GET PHYLO OBJECT
  tree.phylo <- tax_table_to_phylo(~Life/Kingdom/Phylum/Class/Order/Family/Genus,
                                   data = full_tax, collapse = TRUE)
  
  
  ## ***********************
  ## 4 - GET TREEDATA OBJECT
  
  # We will define new labels for the tips
  d <- data.frame(label=full_tax$Genus,
                  new_label=tax[,"Genus"],
                  phylum=tax[,"Phylum"])
  d[!d$phylum %in% c("Acidobacteriota","Actinobacteriota", "Bacteroidota", "Chloroflexi" , "Firmicutes", "Proteobacteria"), "phylum"] <- "Other"
  
  # Create a treedata object, that will contain the new_label
  tree.td <- full_join(tree.phylo, d, by="label")
  
  # Get the list of tip labels belonging to each phylum (phyla are the names of the elements in the list)
  l <- list()
  for (phylum in unique(d$phylum)){
    print(phylum)
    l[[phylum]] <- d[d$phylum==phylum,"label"]
  }
  
  # Group the nodes and edges by phylum (to color the branches by phylum)
  tree.td <- groupOTU(tree.td, l)
  
  
  ## ***********************
  ## 5 - RETURN TREEDATA OBJECT
  return(tree.td)
}


#_________________________________
## 3.3. Function to get dataframe with information on number of ASVs detected for each genus ####
nb_ASV_per_genus <- function(physeq){
  
  asvGenus.df <- as.data.frame(tax_table(physeq)) %>%
    rownames_to_column("ASV") %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    select(ASV, Genus, leaf) %>%
    # get nb of ASV per genus
    group_by(leaf, Genus) %>%
    summarize(n_ASV = n_distinct(ASV)) %>%
    # about 30,851 ASVs belong to an unknown genus
    filter(!is.na(Genus))
  
  asvGenus.df <- data.frame(n_ASV=asvGenus.df$n_ASV,
                            row.names = asvGenus.df$leaf)
  
  # Return df
  return(asvGenus.df)
}


#______________________
# Function to get a dataframe with information on average count per sample for each genus
count_per_genus <- function(physeq){
  # Get relative counts (common-scale normalization: sum per sample==100)
  cat("Common-scale normalization...\n")
  physeq.CSN <- transform_sample_counts(physeq, function(x) (x*100) / sum(x) )
  print(table(sample_sums(physeq.CSN))) # sanity check total of 100 counts per sample
  
  # Get the mean count per sample for each ASV
  cat("\nAverage count of each genus per sample....\n")
  meanCountASV <- t(data.frame(as.list(colMeans(otu_table(physeq.CSN)))))
  # sum(meanCountASV[,1]) # sanity check should sum to 100
  
  # Get a dataframe with the average count of each Genus per sample
  countGenus.df <- meanCountASV %>%
    as.data.frame()%>%
    rownames_to_column("ASV") %>%
    rename(Count=V1) %>%
    # Add taxonomy
    left_join(as.data.frame(tax_table(physeq.CSN)) %>%
                rownames_to_column("ASV"),
              by="ASV") %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    # get average count per sample for each Genus
    group_by(leaf, Genus) %>%
    summarize(count=sum(Count)) %>% # if several ASVs belong to same Genus, need to sum their counts
    # about 30,851 ASVs belong to an unknown genus (which represents ~10% of the counts)
    filter(!is.na(Genus))
  
  countGenus.df <- data.frame(count=countGenus.df$count,
                              row.names = countGenus.df$leaf)
  
  # Return dataframe
  return(countGenus.df)
}


#_________________________________
## 3.4. Function to get a dataframe with information on number of datasets where this genus is present ####
dataset_per_genus <- function(physeq){
  
  print("CAUTION: give the agglomerated phyloseq object for shorter running time")
  # Get a dataframe with the number of dataset(s) each Genus is found in
  datasetGenus.df <- physeq %>%
    # melt to long format
    psmelt() %>%
    # keep only ASVs with a count != 0 in this specific sample
    filter(Abundance != 0) %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    # get the number of datasets with non-0 count for each genus
    group_by(leaf, Genus) %>%
    summarize(n_dataset = n_distinct(Group))
  
  datasetGenus.df <- data.frame(n_dataset=datasetGenus.df$n_dataset,
                                row.names = datasetGenus.df$leaf)
  
  # Return dataframe
  return(datasetGenus.df)
}


# ******************
# 4. PLOT TREES ####
# ******************

# path.data <- file.path(path.root, "data/analysis-combined/01_TaxonomicTree_1")
colors <- paste0(brewer.pal(8, "Dark2"), "99", sep="")
names(colors) <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Firmicutes", "Proteobacteria", "Other")

## *********************
## 4.1. All samples ####

# Agglomerate to Genus level
suppressWarnings(suppressMessages(physeq.glom <- tax_glom(physeq.all_3, "Genus")))

# Get the treedata object, and info on nb of ASV per genus, average count of each genus per sample, sample types containing the genus.
tree.all <- phyloseq_to_treedata(physeq=physeq.glom)

asvDF.all <- nb_ASV_per_genus(physeq=physeq.all_3)

# GroupDF.all <- count_per_genus(physeq=physeq.all)
datasetDF.all <- dataset_per_genus(physeq=physeq.glom)

SubGroup.df <- physeq.glom %>%
  # melt to long format
  psmelt() %>%
  # keep only ASVs with a count != 0
  filter(Abundance != 0) %>%
  # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
  mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
  # get the sample_type(s) each genus can be found in
  group_by(leaf, Genus) %>%
  summarize(spltype = paste(sort(unique(sample_Class)), collapse=","))

SubGroup.df <- data.frame(spltype=SubGroup.df$spltype,
                         row.names=SubGroup.df$leaf)


# Plot tree (all genera)
p1 <- ggtree(tree.all, layout="circular", aes(color=group), lwd=0.4)+
  # geom_tiplab(aes(label=new_label), size=2, color="black")+
  scale_color_manual(values=colors)+
  labs(color="Phylum")


# Add outer layers
p2 <- gheatmap(p1, asvDF.all, offset=0.5, width=.08, colnames=FALSE)+
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = "#2b8cbe",
                       trans="log",
                       breaks=c(1,10,100,1000),
                       name="Number of ASVs")

# p3 <- gheatmap(p2+new_scale_fill(),
#                GroupDF.all, offset=1.1, width=.08, colnames = FALSE)+
#   scale_fill_gradient(low = "white",
#                       high = muted("red"),
#                       trans="log",
#                       breaks=c(1e-5,1e-3,1e-1,10),
#                       labels = parse(text=c("10^-5","10^-3","10^-1", "10")),
#                       name="Mean relative abundance (%)")+
#   labs(title="All samples")

p3 <- gheatmap(p2+new_scale_fill(),
               datasetDF.all, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      breaks = c(1,5,10,14),
                      name="Number of samples")

p4 <- gheatmap(p3 + new_scale_fill(),
               SubGroup.df, offset = 1.7, width = .08, colnames = FALSE) +
  scale_fill_manual(values = c("#66c2a5", "#b3cde3", "#ccebc5", "#decbe4",
                               #"#fed9a6", 
                               "#ffffcc", #"#e5d8bd", "#fddaec",
                               "#fb9a99", "#b2df8a")
                               , #"#999999"),  # added 3 colors
                    name = "Samples") +
  labs(title = "(a) Phylogenetic tree of soil microbiome") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0)) +
  scale_y_reverse()
p4

ggsave(p4, file = "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Taxonomy/taxTree_all_without_seed_control_soil.svg", width=12, height=8, dpi = 1000)


##Plot Specific phylum
library(phyloseq)
library(microbiomeutilities)
library(ggpubr)
library(patchwork)
library(scales)

# Load phyloseq object
suppressWarnings(suppressMessages(ps_glom <- tax_glom(physeq.all_3, "Genus")))
suppressWarnings(suppressMessages(p0.f <- format_to_besthit(ps_glom)))

# # Selected ASVs
# taxa_list <- c("36fd004179a84695b3e40f1ce2490cd8:c__Actinobacteriota",
#                "37cabbdb27757ee34a11939985c7f668:c__Proteobacteria",
#                "099c69de18d95328168346332ef4288b:c__Firmicutes",
#                "04f3a0aaca3b674b39d533131f1ac888:c__Bacteroidota"
# )

# Panel labels
panel_labels <- paste0("(", letters[1:10], ".)")

# Group order & colors
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
  # "Seed"
)



# mycols <- c("brown3", "steelblue", "grey50")

mycols <- c("#E69F00",  # orange
            "#56B4E9",  # sky blue
            "#009E73",  # bluish green
            "#F0E442",  # yellow
            "#0072B2",  # blue
            "#D55E00",  # vermilion
            "#CC79A7",  # reddish purple
            "#999999",  # gray
            "#e41a1c",  # red
            "#4daf4a",  # green
            "#984ea3",  # purple
            "#ff7f00",  # orange
            "#a65628",  # brown
            "#f781bf") # pink

pn <- plot_taxa_boxplot(p0.f,
                        taxonomic.level = "Phylum",
                        top.otu = 6, 
                        group = "Group",
                        add.violin= FALSE,
                        title = "Top four phylum", 
                        keep.other = FALSE,
                        group.order = group_order,
                        group.colors = mycols,
                        dot.size = 1) + 
  theme(
    axis.text  = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 8, face = "bold"),
    axis.line  = element_line(colour = "black"),
    panel.background = element_blank(),
    legend.text = element_text(size = 8),
    legend.direction="vertical",
    legend.position="none"
  )

pn
ggsave(pn, file = "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Taxonomy/combined_plot.svg", height = 8, width = 14, dpi = 1000)

library(MicrobiotaProcess)
mpse4 <- ps_glom %>% as.MPSE() 

# visualize the relative abundance of top 20 phyla for each sample.
p1 <- mpse4 %>%
  mp_plot_abundance(
    .abundance=Abundance,
    .group=SubGroup, 
    taxa.class = Genus, 
    topn = 20,
    relative = TRUE
  )

mpse4$SubGroup <- factor(
  mpse4$SubGroup,
  # levels = c("Control_Soil", "Bulk_Soil", "Rhizosphere_Soil")
  levels = c("Bulk_Soil", "Rhizosphere_Soil")
)

p1 <- mpse4 %>%
  mp_plot_abundance(
    .abundance = Abundance,
    .group = Class,
    taxa.class = Genus,
    topn = 20,
    relative = TRUE,
    plot.group = F
  )
p1
ggsave(p1, file = "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Taxonomy/amplicon_top_20.svg", height = 6, width = 14)



mpse4$Group <- factor(
  mpse4$Group,
  levels = c(
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
             "21D_AB211RE")
)

p2 <- mpse4 %>%
  mp_plot_abundance(
    .abundance = Abundance,
    .group = Group,
    taxa.class = Genus,
    topn = 20,
    relative = TRUE,
    plot.group = T,
    
  )
p2
ggsave(p2, file = "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/Final_Result/Amplicon/Taxonomy/amplicon_top_20_group.svg", height = 5, width = 8)

