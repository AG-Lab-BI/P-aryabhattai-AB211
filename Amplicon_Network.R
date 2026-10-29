# Load necessary libraries
library(file2meco)
library(microeco)
library(phyloseq)

set.seed(124)
# Set working directory (update the path accordingly)
setwd("/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/qiime2_biom_files")

# Input files
otuqzafile <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/qiime2_biom_files/tableNoFilt.qza"
taxaqzafile <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/qiime2_biom_files/taxonomyNoFilt.qza"
mapfile <- "/home/sslab4/Documents/Dibakar_Roy/Collaboration/Sangita_Abhrajyoti_Sir/qiime2_biom_files/metadata.tsv"
tree_data <- "treeNoFilt.qza"
rep_data <- "repseqsNoFilt.qza"

test1 <- qiime2meco(otuqzafile, sample_table = mapfile, taxonomy_table = taxaqzafile, phylo_tree = tree_data, rep_fasta = rep_data, auto_tidy = TRUE)
test1
# ps_tree@tax_table@.Data = substring(ps_tree@tax_table@.Data, 4)

suppressWarnings(suppressMessages(ps <- meco2phyloseq(test1)))
ps@tax_table@.Data = substring(ps@tax_table@.Data, 4)


physeq.all <- prune_taxa(taxa_sums(ps)>0, ps) # remove ASVs that are not present anymore
physeq.all <- subset_taxa(physeq.all, Kingdom != "Eukaryota")
physeq.all_1 <- subset_taxa(physeq.all, Family != "Mitochondria")
physeq.all_2 <- subset_taxa(physeq.all_1, Order != "Chloroplast")

ps_glom <- tax_glom(physeq.all_2, "Genus", NArm = TRUE)

cat("Nb of Sputum samples:", nsamples(ps))


D0_Control <- subset_samples(ps_glom, SubGroup == "0 Day")

meco_D0_Control<- phyloseq2meco(D0_Control)
meco_D0_Control$tidy_dataset()

# Create the network
t1 <- trans_network$new(
  dataset = meco_D0_Control,
  taxa_level = "Genus",
  cor_method = "sparcc", 
  use_sparcc_method = "SpiecEasi", 
  filter_thres = 0.0001, 
  nThreads = 16
)

# Calculate the network
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Perform clustering (community detection)
t1$cal_module(method = "cluster_fast_greedy")

# Save the network to a GEXF file
t1$save_network(filepath = "0 Day.gexf")



D0_AB211 <- subset_samples(ps_glom, SubGroup == "14 Day")

meco_D0_AB211 <- phyloseq2meco(D0_AB211)
meco_D0_AB211$tidy_dataset()

# Create the network
t2 <- trans_network$new(
  dataset = meco_D0_AB211,
  taxa_level = "Genus",
  cor_method = "sparcc", 
  use_sparcc_method = "SpiecEasi", 
  filter_thres = 0.0001, 
  nThreads = 12
)

# Calculate the network
t2$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Perform clustering (community detection)
t2$cal_module(method = "cluster_fast_greedy")

# Save the network to a GEXF file
t2$save_network(filepath = "14 Day.gexf")



D0_AB211_RE <- subset_samples(ps_glom, SubGroup == "7 Day")

meco_D0_AB211_RE <- phyloseq2meco(D0_AB211_RE)
meco_D0_AB211_RE$tidy_dataset()

# Create the network
t3 <- trans_network$new(
  dataset = meco_D0_AB211_RE,
  taxa_level = "Genus",
  cor_method = "sparcc", 
  use_sparcc_method = "SpiecEasi", 
  filter_thres = 0.0001, 
  nThreads = 12
)

# Calculate the network
t3$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Perform clustering (community detection)
t3$cal_module(method = "cluster_fast_greedy")

# Save the network to a GEXF file
t3$save_network(filepath = "7 Day.gexf")




D7_Control <- subset_samples(ps_glom, SubGroup == "21 Day")

meco_D7_Control <- phyloseq2meco(D7_Control)
meco_D7_Control$tidy_dataset()

# Create the network
t4 <- trans_network$new(
  dataset = meco_D7_Control,
  taxa_level = "Genus",
  cor_method = "sparcc", 
  use_sparcc_method = "SpiecEasi", 
  filter_thres = 0.0001, 
  nThreads = 12
)

# Calculate the network
t4$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Perform clustering (community detection)
t4$cal_module(method = "cluster_fast_greedy")

# Save the network to a GEXF file
t4$save_network(filepath = "21 Day.gexf")


