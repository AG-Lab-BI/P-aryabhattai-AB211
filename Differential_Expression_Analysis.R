#################################################
DESeq2
########################
library(readxl)
#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(DESeq2)
library(writexl)
library(tibble)
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/trimmed_count_matrix.xlsx"
setwd(dirname(PATH))
##read the data as dataframe and first column gene
DATA = data.frame(read_excel(PATH, sheet = "OD1"), row.names = 1)
view(DATA)
##remove genes with 0 values or
DATA = DATA[rowSums(DATA !=0) > 0, ]
META = data.frame(read_excel(PATH, sheet = "META_Data_OD1"), row.names = 1)
print(META)
rownames(META)
##rownames of META Data as column names of DATA
MATRIX = DATA[ ,rownames(META)]
view(MATRIX)
dds = DESeqDataSetFromMatrix(countData = MATRIX,
                             colData = META,
                             design = ~Group)
##reorder my groups
dds$Group = relevel(dds$Group, ref = "Control")
dds = DESeq(dds)
RES_DF = data.frame(results(dds))
RES_DF$pvalue[RES_DF$pvalue == 0] <- 1e-300
RES_DF$padj[RES_DF$padj == 0] <- 1e-300
##convert rowname as column using tibble package
RES_DF = rownames_to_column(RES_DF, var = "Gene")
RES_DF[is.na(RES_DF)] = 0
SIG = subset(RES_DF, RES_DF$pvalue <= 0.05)
UP_Gene = subset(SIG, SIG$log2FoldChange >= 1)
DOWN_Gene = subset(SIG, SIG$log2FoldChange <= -1)
RES_LIST = list(ALL_DEG = RES_DF, significant = SIG, Upregulated = UP_Gene, Downregulated = DOWN_Gene)
write_xlsx(RES_LIST, "OD1_LFC1_26032026.xlsx")
write_xlsx(RES_LIST, "OD1_LFC1_26032026.xlsx")
library(readxl)
library(ggplot2)
library(tibble)
library(ggrepel)
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD1_LFC1_26032026.xlsx"
setwd(dirname(PATH))
list.files()
DATA = data.frame(read_excel(PATH, sheet = "ALL_DEG"))
TITLE = "OD1"
DF = DATA[, c("Gene", "log2FoldChange", "pvalue")]
colnames(DF) = c("Gene", "log2FC", "Pvalue")
DF$Expression = "Not-Significant"
DF$Expression[DF$log2FC >= 1 & DF$Pvalue <= 0.05] = "Up-Regulated"
DF$Expression[DF$log2FC <= -1 & DF$Pvalue <= 0.05] = "Down-Regulated"
UP = subset(DF, DF$Expression == "Up-Regulated")
DOWN = subset(DF, DF$Expression == "Down-Regulated")
TOP_UP = head(UP[order(UP$Pvalue), ], 15)
TOP_DOWN = head(DOWN[order(DOWN$Pvalue), ], 15)
LABEL_TO_SHOW = c(TOP_UP$Gene, TOP_DOWN$Gene)
DF$Label = ifelse(DF$Gene %in% LABEL_TO_SHOW, DF$Gene, NA)
PLOT = ggplot(DF, aes(x = log2FC, y = -log10(Pvalue), col = Expression, label = Label))+
  geom_point(size = 1.5)+
  scale_color_manual(values = c("darkblue", "gray", "coral"))+
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.7, linetype = "dashed")+
  geom_vline(xintercept = c(-1.5,1.5), color = "black", linewidth = 0.7, linetype = "dashed")+
  labs(title = TITLE, x = "log2FoldChange", y = "-log10(P-Value)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 10, color = "black", face = "bold"),
        legend.title  = element_text(size = 12, color = "black", face = "bold"))+
  geom_text_repel(size = 3, col = "black", fontface = "bold",
                  min.segment.length = 0.1,
                  point.size = 2,
                  segment.size = 0.5,
                  max.overlaps = Inf)+
  coord_cartesian(ylim = c(0, 400), xlim = c(-8, 8))
PLOT
svg(paste0(TITLE, "OD1_volcano.svg"), height = 6, width = 8)
print(PLOT)
dev.off()
##read the data as dataframe and first column gene
DATA_2 = data.frame(read_excel(PATH, sheet = "OD2"), row.names = 1)
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/trimmed_count_matrix.xlsx"
setwd(dirname(PATH))
##read the data as dataframe and first column gene
DATA_2 = data.frame(read_excel(PATH, sheet = "OD2"), row.names = 1)
view(DATA_2)
##remove genes with 0 values or
DATA_2 = DATA_2[rowSums(DATA !=0) > 0, ]
META_2 = data.frame(read_excel(PATH, sheet = "META_Data_OD2"), row.names = 1)
print(META_2)
rownames(META_2)
##rownames of META Data as column names of DATA
MATRIX_2 = DATA_2[ ,rownames(META_2)]
view(MATRIX_2)
dds = DESeqDataSetFromMatrix(countData = MATRIX_2,
                             colData = META_2,
                             design = ~Group)
##reorder my groups
dds$Group = relevel(dds$Group, ref = "Control")
dds = DESeq(dds)
RES_DF = data.frame(results(dds))
RES_DF$pvalue[RES_DF$pvalue == 0] <- 1e-300
RES_DF$padj[RES_DF$padj == 0] <- 1e-300
##convert rowname as column using tibble package
RES_DF = rownames_to_column(RES_DF, var = "Gene")
RES_DF[is.na(RES_DF)] = 0
SIG = subset(RES_DF, RES_DF$pvalue <= 0.05)
UP_Gene = subset(SIG, SIG$log2FoldChange >= 1)
DOWN_Gene = subset(SIG, SIG$log2FoldChange <= -1)
RES_LIST = list(ALL_DEG = RES_DF, significant = SIG, Upregulated = UP_Gene, Downregulated = DOWN_Gene)
write_xlsx(RES_LIST, "OD2_LFC1_26032026.xlsx")
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD2_LFC1_26032026.xlsx"
setwd(dirname(PATH))
list.files()
DATA = data.frame(read_excel(PATH, sheet = "ALL_DEG"))
TITLE = "OD2"
DF = DATA[, c("Gene", "log2FoldChange", "pvalue")]
colnames(DF) = c("Gene", "log2FC", "Pvalue")
DF$Expression = "Not-Significant"
DF$Expression[DF$log2FC >= 1 & DF$Pvalue <= 0.05] = "Up-Regulated"
DF$Expression[DF$log2FC <= -1 & DF$Pvalue <= 0.05] = "Down-Regulated"
UP = subset(DF, DF$Expression == "Up-Regulated")
DOWN = subset(DF, DF$Expression == "Down-Regulated")
TOP_UP = head(UP[order(UP$Pvalue), ], 15)
TOP_DOWN = head(DOWN[order(DOWN$Pvalue), ], 15)
LABEL_TO_SHOW = c(TOP_UP$Gene, TOP_DOWN$Gene)
DF$Label = ifelse(DF$Gene %in% LABEL_TO_SHOW, DF$Gene, NA)
PLOT = ggplot(DF, aes(x = log2FC, y = -log10(Pvalue), col = Expression, label = Label))+
  geom_point(size = 1.5)+
  scale_color_manual(values = c("darkblue", "gray", "coral"))+
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.7, linetype = "dashed")+
  geom_vline(xintercept = c(-1.5,1.5), color = "black", linewidth = 0.7, linetype = "dashed")+
  labs(title = TITLE, x = "log2FoldChange", y = "-log10(P-Value)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 10, color = "black", face = "bold"),
        legend.title  = element_text(size = 12, color = "black", face = "bold"))+
  geom_text_repel(size = 3, col = "black", fontface = "bold",
                  min.segment.length = 0.1,
                  point.size = 2,
                  segment.size = 0.5,
                  max.overlaps = Inf)+
  coord_cartesian(ylim = c(0, 400), xlim = c(-8, 8))
PLOT
svg(paste0(TITLE, "OD2_volcano.svg"), height = 6, width = 8)
print(PLOT)
dev.off()
setwd("Results")
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/trimmed_count_matrix.xlsx"
OD1 = read_excel("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD1_LFC1_26032026.xlsx", sheet = 1)
OD2 = read_excel("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD2_LFC1_26032026.xlsx", sheet = 1)
MERGED_DF  = Reduce(function(x,y) merge(x, y, by = "Gene", all = TRUE), list(OD1, OD2))
view(MERGED_DF)
colnames(MERGED_DF)
MERGED_DF$Avg_Pval = (MERGED_DF$pvalue.x + MERGED_DF$pvalue.y)/2
colnames(MERGED_DF)
PVAL_COL = "pvalue_OD1"
logFC_COL = "log2FoldChange_OD1"
SIG_MERGED = subset(MERGED_DF,MERGED_DF[PVAL_COL] <= 0.05 )
colnames(MERGED_DF)
PVAL_COL = "pvalue.x"
logFC_COL = "log2FoldChange.x"
SIG_MERGED = subset(MERGED_DF,MERGED_DF[PVAL_COL] <= 0.05 )
UP_MERGED  = subset(SIG_MERGED, SIG_MERGED[logFC_COL] >= 1)
DOWN_MERGED  = subset(SIG_MERGED, SIG_MERGED[logFC_COL] <= -1)
view(UP_MERGED)
ALL_RES_LIST  = list(ALL_DEG = MERGED_DF, Significant = SIG_MERGED, UpRegulated = UP_MERGED, DownRegulated = DOWN_MERGED)
write_xlsx(ALL_RES_LIST, "Combined_DEG_Result.xlsx")
view(MERGED_DF)
setwd("Results")
PATH = "E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/trimmed_count_matrix.xlsx"
OD1 = read_excel("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD1_LFC1_26032026.xlsx", sheet = 1)
OD2 = read_excel("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/OD2_LFC1_26032026.xlsx", sheet = 1)
MERGED_DF  = Reduce(function(x,y) merge(x, y, by = "Gene", all = TRUE), list(OD1, OD2))
view(MERGED_DF)
write_xlsx(MERGED_DF, "Merged.xlsx")

##########################################
####Heatmap plots
###################
DATA = data.frame(read_xlsx("AB211_KO_sorted.xlsx", sheet = "Siderophore_LFC_1"),row.names = 1)
colnames(DATA)
MATRIX = DATA[, c("log2FoldChange_OD1", "log2FoldChange_OD2")]
colnames(DATA) = c("Exponential phase", "Stationary phase")
PLOT_DATA = Heatmap(as.matrix(DATA),
                    name = "log2FC",
                    col = colorRamp2(c(-2, 0, 2), c("darkblue","yellow", "darkred")),
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    row_names_gp = gpar(fontface = "bold.italic",fontfamily = "Arial"),
                    column_names_gp = gpar(fontface = "bold",fontfamily = "Arial"))
PLOT_DATA
svglite("Siderophore_LFC_1.svg", height = 8, width =8)
##png("Heatmap_amino.png", height = 8, width =8, units = "in", res = 300)
print(PLOT_DATA)
dev.off()
DATA = data.frame(read_xlsx("AB211_KO_sorted.xlsx", sheet = "Flagella_motility_LFC_1"),row.names = 1)
colnames(DATA)
MATRIX = DATA[, c("log2FoldChange_OD1", "log2FoldChange_OD2")]
colnames(DATA) = c("Exponential phase", "Stationary phase")
PLOT_DATA = Heatmap(as.matrix(DATA),
                    name = "log2FC",
                    col = colorRamp2(c(-2, 0, 2), c("darkblue","yellow", "darkred")),
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    row_names_gp = gpar(fontface = "bold.italic",fontfamily = "Arial"),
                    column_names_gp = gpar(fontface = "bold",fontfamily = "Arial"))
PLOT_DATA
svglite("Flagella_motility_LFC_1.svg", height = 8, width =8)
##png("Heatmap_amino.png", height = 8, width =8, units = "in", res = 300)
print(PLOT_DATA)
dev.off()



###############################################################
##KEGG_enrichment
###################################################################
## R commands
###################################################################
library(readxl)
###################################################################
## Set working directory
###################################################################
setwd("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/Results/")
###################################################################
## Taking input file and writing a list file
###################################################################
final_input_file<-read_xlsx("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/Results/KEGG_enrichment_list.xlsx", sheet = "Exponential_LFC1_0.05") # WITH HEADER AS "ID" & "FC"
write.table(final_input_file$ID, quote = F, row.names = F, col.names = F, "list_OD1.txt")


###################################################################
## BASH commands
###################################################################
# Known = Priestia megaterium QM B1551
# Unknown = Priestia aryabhatti AB211
# ###################################################################
# # Download known gene seq in translated cds (locus tag) B1551 (Genbank only for ID compatibility issues)
# ## bash command
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/025/825/GCA_000025825.1_ASM2582v1/GCA_000025825.1_ASM2582v1_translated_cds.faa.gz
# # Download unknown gene seq in translated cds (locus tag) AB211
# # bash command
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/395/GCF_001858395.1_ASM185839v1/GCF_001858395.1_ASM185839v1_translated_cds.faa.gz
###################################################################
# Extract the files
###################################################################
## bash command
# for f in *.gz; do gunzip $f; done
####################################################################
# Shorten FASTA header
# sed -n '/^>/ {s/.*locus_tag=\([^]]*\).*/>\1/; p; n}; p' GCA_000025825.1_ASM2582v1_translated_cds.faa > B1551.fasta
# sed -n '/^>/ {s/.*locus_tag=\([^]]*\).*/>\1/; p; n}; p' GCF_001858395.1_ASM185839v1_translated_cds.faa > AB211.fasta
####################################################################
# Add 'RS' after hyphen in list file
# sed 's/_\([0-9]\+\)/_RS\1/' list_OD1.txt > modified_list_OD1.txt


###################################################################
## R commands
###################################################################
# Take the modified input in a dataframe
final_input_modified<-read.table("list_OD1.txt") # WITH HEADER AS "ID" & "FC"
colnames(final_input_modified)[1] <- c("AB211")
## Extract only the sequence from AB211 of ids given in the modified_list.txt file

# Load necessary library
library(Biostrings)
# Step 1: Read the list of IDs
ids <- readLines("list_OD1.txt")
# Step 2: Read the FASTA file
fasta_sequences <- readAAStringSet("AB211_GCA.fasta")
# Step 3: Extract the matched sequences
matched_sequences <- fasta_sequences[names(fasta_sequences) %in% ids]
# Step 4: Write the matched sequences to a new FASTA file
writeXStringSet(matched_sequences, filepath = "AB211_OD1.fasta")

###################################################################
## BASH commands
###################################################################
## Do BLAST
# Creating BLAST database (BASH)
# makeblastdb -in AB211_modified_OD1.fasta -dbtype prot -out known
# blastp -query B1551.fasta -db known -evalue 1e-5 -outfmt 6 -out results_ab211_b1551_OD1.txt
## we got 248 hits in 1e-5 and 266 hits in 1e-5
###################################################################

###################################################################
## R commands
###################################################################
# Input blast result
blast_input_file<-read.table("results_ab211_b1551_OD1.txt", header =F)
# remove duplicates hits from database entry
blast_input_file_unique <- blast_input_file[!duplicated(blast_input_file[, 2]), ]
# Keep only the first two columns
blast_IDs <- blast_input_file_unique[, 1:2]
# Change the headers
colnames(blast_IDs)[1:2] <- c("B1551", "AB211")
# Add "RS" after every underscore in the AB211 column
#final_input_file$ID <- gsub("_", "_RS", final_input_file$ID)
# Perform the merge operation to match IDs and get the corresponding FC values
matched_data <- merge(blast_IDs, final_input_file, by.x = "AB211", by.y = "ID")
# Select only the matched ID and FC values
final_output <- matched_data[, c("B1551", "FC")]
colnames(final_output)[1] <- c("ID")
# Remove "RS" after every underscore in the AB211 column
#final_output$ID <- gsub("_RS", "_", final_output$ID)
#########################################################################

BiocManager::install("clusterProfiler")


library(clusterProfiler)
#########################################################################
# GSEA analysis requires a ranked gene list, which contains three features:
# numeric vector: fold change or other type of numerical variable
# named vector: every number has a name, the corresponding gene ID
# sorted vector: number should be sorted in decreasing order

## feature 1: numeric vector
geneList <- final_output[,2]
## feature 2: named vector
names(geneList) <- as.character(final_output[,1])
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)

## KEGG pathway over-representation analysis
# see https://rest.kegg.jp/list/organism
kk <- enrichKEGG(gene         = gene,
                 organism     = 'bmq',
                 pvalueCutoff = 0.05)
dotplot(kk, showCategory = kk$Description[1:10])
## No enriched KEGG category, we need to show it differently
final_KEGG_result<-kk@result
# Remove repetitive word
final_KEGG_result$Description <- gsub("- Priestia megaterium QM B1551", "", final_KEGG_result$Description)

# Split the values into numerator and denominator and compute the ratio for Gene
final_KEGG_result$Gratio <- sapply(strsplit(as.character(final_KEGG_result$GeneRatio), "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})
# Split the values into numerator and denominator and compute the ratio for Background
final_KEGG_result$Bratio <- sapply(strsplit(as.character(final_KEGG_result$BgRatio), "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})
library(writexl)
write_xlsx(final_KEGG_result, "AB211_bmq_KEGG_OD1.xlsx")


## ID CONVERSION (Locus tag to protein ID)
###################################################################
## BASH commands
###################################################################
## Update ID according to current Refseq notation
final_KEGG_result$ID <- gsub("bmq", "BMQ_RS", final_KEGG_result$ID)
## Download the refseq file with both IDs (locus tag and protein)
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/825/GCF_000025825.1_ASM2582v1/GCF_000025825.1_ASM2582v1_translated_cds.faa.gz
## Extract
# for f in *.gz; do gunzip $f; done

###################################################################
## R commands
###################################################################
# Load necessary library
#if(!require(stringr)) install.packages("stringr", dependencies=TRUE)
library(stringr)

# Read the input FASTA file
input_file <- "GCA_000025825.1_ASM2582v1_translated_cds.faa"
output_file <- "output.txt"

# Read the file lines
lines <- readLines(input_file)

# Initialize vectors to store locus_tag and protein_id
locus_tags <- c()
protein_ids <- c()

# Loop through the lines and extract locus_tag and protein_id
for (line in lines) {
  if (startsWith(line, ">")) {
    # Extract locus_tag
    locus_tag <- str_extract(line, "\\[locus_tag=([^]]+)\\]")
    locus_tag <- str_replace_all(locus_tag, "\\[locus_tag=|\\]", "")
    locus_tags <- c(locus_tags, locus_tag)
    
    # Extract protein_id
    protein_id <- str_extract(line, "\\[protein_id=([^]]+)\\]")
    protein_id <- str_replace_all(protein_id, "\\[protein_id=|\\]", "")
    protein_ids <- c(protein_ids, protein_id)
  }
}
library(readxl)
library(dplyr)
library(tidyr)

kegg <- read_excel("AB211_bmq_KEGG_OD1.xlsx")
kegg_long <- kegg %>%
  separate_rows(geneID, sep = "/")
# Create a data frame
locus_protein_IDs <- data.frame(Locus_tag = locus_tags, Protein_ID = protein_ids)

# Get the matching IDs
final_KEGG_result_pIDs <- merge(kegg_long, locus_protein_IDs, by.x = "geneID", by.y = "Locus_tag")
## Notice that this has several less IDS than the original final_KEGG_result,
## that's because the ommitted ones are the RNA's (t-RNA, r-RNA) that do have
## locus tag but no proteins

# FC of Gene to Background Ratio
final_KEGG_result_pIDs$FC <- final_KEGG_result_pIDs$Gratio / final_KEGG_result_pIDs$Bratio
#################
## Plot by p-value
#################
# Sort the data frame by the 'pvalue' column in ascending order
library(dplyr)
final_KEGG_result_sorted <- final_KEGG_result_pIDs %>% arrange(pvalue)
write_xlsx(final_KEGG_result_sorted, "sorted_OD1.xlsx")
sorted_OD1 <- read_excel("sorted_OD1.xlsx")
# Keep only the first 12 rows
final_KEGG_result_sorted_first12 <- head(sorted_OD1, 12)

# Ensure FC and pvalue are numeric
final_KEGG_result_sorted_first12$FC <- as.numeric(final_KEGG_result_sorted_first12$FC)
final_KEGG_result_sorted_first12$pvalue <- as.numeric(final_KEGG_result_sorted_first12$pvalue)





##########################
###OD2 KEGG Enrichment


###################################################################
## R commands
###################################################################
library(readxl)
###################################################################
## Set working directory
###################################################################
setwd("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/Results/")
###################################################################
## Taking input file and writing a list file
###################################################################
final_input_file<-read_xlsx("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/Results/KEGG_enrichment_list.xlsx", sheet = "Stationary_LFC1_0.05") # WITH HEADER AS "ID" & "FC"
write.table(final_input_file$ID, quote = F, row.names = F, col.names = F, "list_OD2.txt")


###################################################################
## BASH commands
###################################################################
# Known = Priestia megaterium QM B1551
# Unknown = Priestia aryabhatti AB211
# ###################################################################
# # Download known gene seq in translated cds (locus tag) B1551 (Genbank only for ID compatibility issues)
# ## bash command
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/025/825/GCA_000025825.1_ASM2582v1/GCA_000025825.1_ASM2582v1_translated_cds.faa.gz
# # Download unknown gene seq in translated cds (locus tag) AB211
# # bash command
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/395/GCF_001858395.1_ASM185839v1/GCF_001858395.1_ASM185839v1_translated_cds.faa.gz
###################################################################
# Extract the files
###################################################################
## bash command
# for f in *.gz; do gunzip $f; done
####################################################################
# Shorten FASTA header
# sed -n '/^>/ {s/.*locus_tag=\([^]]*\).*/>\1/; p; n}; p' GCA_000025825.1_ASM2582v1_translated_cds.faa > B1551.fasta
# sed -n '/^>/ {s/.*locus_tag=\([^]]*\).*/>\1/; p; n}; p' GCF_001858395.1_ASM185839v1_translated_cds.faa > AB211.fasta
####################################################################
# Add 'RS' after hyphen in list file
# sed 's/_\([0-9]\+\)/_RS\1/' list_OD1.txt > modified_list_OD1.txt


###################################################################
## R commands
###################################################################
# Take the modified input in a dataframe
final_input_modified<-read.table("list_OD2.txt") # WITH HEADER AS "ID" & "FC"
colnames(final_input_modified)[1] <- c("AB211")
## Extract only the sequence from AB211 of ids given in the modified_list.txt file

# Load necessary library
library(Biostrings)
# Step 1: Read the list of IDs
ids <- readLines("list_OD2.txt")
# Step 2: Read the FASTA file
fasta_sequences <- readAAStringSet("AB211_GCA.fasta")
# Step 3: Extract the matched sequences
matched_sequences <- fasta_sequences[names(fasta_sequences) %in% ids]
# Step 4: Write the matched sequences to a new FASTA file
writeXStringSet(matched_sequences, filepath = "AB211_OD2.fasta")

###################################################################
## BASH commands
###################################################################
## Do BLAST
# Creating BLAST database (BASH)
# makeblastdb -in AB211_OD2.fasta -dbtype prot -out known
# blastp -query B1551.fasta -db known -evalue 1e-5 -outfmt 6 -out results_ab211_b1551_OD2.txt
## we got 248 hits in 1e-5 and 266 hits in 1e-5
###################################################################

###################################################################
## R commands
###################################################################
# Input blast result
blast_input_file<-read.table("results_ab211_b1551_OD2.txt", header =F)
# remove duplicates hits from database entry
blast_input_file_unique <- blast_input_file[!duplicated(blast_input_file[, 2]), ]
# Keep only the first two columns
blast_IDs <- blast_input_file_unique[, 1:2]
# Change the headers
colnames(blast_IDs)[1:2] <- c("B1551", "AB211")
# Add "RS" after every underscore in the AB211 column
#final_input_file$ID <- gsub("_", "_RS", final_input_file$ID)
# Perform the merge operation to match IDs and get the corresponding FC values
matched_data <- merge(blast_IDs, final_input_file, by.x = "AB211", by.y = "ID")
# Select only the matched ID and FC values
final_output <- matched_data[, c("B1551", "FC")]
colnames(final_output)[1] <- c("ID")
# Remove "RS" after every underscore in the AB211 column
#final_output$ID <- gsub("_RS", "_", final_output$ID)
#########################################################################

BiocManager::install("clusterProfiler")


library(clusterProfiler)
#########################################################################
# GSEA analysis requires a ranked gene list, which contains three features:
# numeric vector: fold change or other type of numerical variable
# named vector: every number has a name, the corresponding gene ID
# sorted vector: number should be sorted in decreasing order

## feature 1: numeric vector
geneList <- final_output[,2]
## feature 2: named vector
names(geneList) <- as.character(final_output[,1])
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)

## KEGG pathway over-representation analysis
# see https://rest.kegg.jp/list/organism
kk <- enrichKEGG(gene         = gene,
                 organism     = 'bmq',
                 pvalueCutoff = 0.05)
dotplot(kk, showCategory = kk$Description[1:10])
## No enriched KEGG category, we need to show it differently
final_KEGG_result<-kk@result
# Remove repetitive word
final_KEGG_result$Description <- gsub("- Priestia megaterium QM B1551", "", final_KEGG_result$Description)

# Split the values into numerator and denominator and compute the ratio for Gene
final_KEGG_result$Gratio <- sapply(strsplit(as.character(final_KEGG_result$GeneRatio), "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})
# Split the values into numerator and denominator and compute the ratio for Background
final_KEGG_result$Bratio <- sapply(strsplit(as.character(final_KEGG_result$BgRatio), "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})
library(writexl)
write_xlsx(final_KEGG_result, "AB211_bmq_KEGG_OD2.xlsx")


## ID CONVERSION (Locus tag to protein ID)
###################################################################
## BASH commands
###################################################################
## Update ID according to current Refseq notation
#final_KEGG_result$ID <- gsub("bmq", "BMQ_RS", final_KEGG_result$ID)
## Download the refseq file with both IDs (locus tag and protein)
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/825/GCF_000025825.1_ASM2582v1/GCF_000025825.1_ASM2582v1_translated_cds.faa.gz
## Extract
# for f in *.gz; do gunzip $f; done

###################################################################
## R commands
###################################################################
# Load necessary library
#if(!require(stringr)) install.packages("stringr", dependencies=TRUE)
library(stringr)

# Read the input FASTA file
input_file <- "GCA_000025825.1_ASM2582v1_translated_cds.faa"
output_file <- "output.txt"

# Read the file lines
lines <- readLines(input_file)

# Initialize vectors to store locus_tag and protein_id
locus_tags <- c()
protein_ids <- c()

# Loop through the lines and extract locus_tag and protein_id
for (line in lines) {
  if (startsWith(line, ">")) {
    # Extract locus_tag
    locus_tag <- str_extract(line, "\\[locus_tag=([^]]+)\\]")
    locus_tag <- str_replace_all(locus_tag, "\\[locus_tag=|\\]", "")
    locus_tags <- c(locus_tags, locus_tag)
    
    # Extract protein_id
    protein_id <- str_extract(line, "\\[protein_id=([^]]+)\\]")
    protein_id <- str_replace_all(protein_id, "\\[protein_id=|\\]", "")
    protein_ids <- c(protein_ids, protein_id)
  }
}
library(readxl)
library(dplyr)
library(tidyr)

kegg <- read_excel("AB211_bmq_KEGG_OD2.xlsx")
kegg_long <- kegg %>%
  separate_rows(geneID, sep = "/")
# Create a data frame
locus_protein_IDs <- data.frame(Locus_tag = locus_tags, Protein_ID = protein_ids)

# Get the matching IDs
final_KEGG_result_pIDs <- merge(kegg_long, locus_protein_IDs, by.x = "geneID", by.y = "Locus_tag")
## Notice that this has several less IDS than the original final_KEGG_result,
## that's because the ommitted ones are the RNA's (t-RNA, r-RNA) that do have
## locus tag but no proteins

# FC of Gene to Background Ratio
final_KEGG_result_pIDs$FC <- final_KEGG_result_pIDs$Gratio / final_KEGG_result_pIDs$Bratio
#################
## Plot by p-value
#################
# Sort the data frame by the 'pvalue' column in ascending order
library(dplyr)
final_KEGG_result_sorted <- final_KEGG_result_pIDs %>% arrange(pvalue)
write_xlsx(final_KEGG_result_sorted, "sorted_OD2.xlsx")
#Open sorted_OD2.xlsx in excel and remove duplicates from the "ID" column i.e bmq0001 column and save it as sorted_OD2.xlsx
sorted_OD2 <- read_excel("sorted_OD2.xlsx")
# Keep only the first 12 rows
final_KEGG_result_sorted_first12 <- head(sorted_OD2, 12)

# Ensure FC and pvalue are numeric
final_KEGG_result_sorted_first12$FC <- as.numeric(final_KEGG_result_sorted_first12$FC)
final_KEGG_result_sorted_first12$pvalue <- as.numeric(final_KEGG_result_sorted_first12$pvalue)

###### GSEA bar plot ######
#OD1
library(readxl)
library(ggplot2)
library(dplyr)
setwd("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/")
pathway_data <- read_xlsx("E:/Ssangita_backup/A_Ghosh/2026/KEGG_enrichment_list.xlsx", sheet ="Sorted_OD1")

# Check the data
head(pathway_data)


# Convert Description to factor & order by Count
data <- pathway_data %>% arrange(desc(Count))

# Plot the data with both p-value (color) and FC (bar border)
ggplot(data, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_bar(stat = "identity", color = "black", aes(size = abs(FC))) +  # Bar outline thickness based on FC
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red", name = "P-value") +
  scale_size_continuous(name = "Fold Change", range = c(0.2, 2)) +  # Controls thickness
  labs(title = "Pathway Enrichment Analysis",
       x = "Pathway",
       y = "Gene Count") +
  theme_minimal(base_size = 14)

library(ggplot2)
library(dplyr)

# Sort by Gene Count for better visual representation
df <- data %>%
  arrange(desc(Count))
library(dplyr)
sorted_FC <- data %>% arrange(desc(FC))
# Filter the data frame to include only rows with 'FC' values >= 1.5
final_KEGG_result_sorted_FC <- sorted_FC %>% filter(FC >= 1.5)
# Plotting
OD1 <- ggplot(final_KEGG_result_sorted_FC, aes(x=reorder(Description, Count), y=Count, fill=pvalue)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low="magenta4", high="coral", name="P-value") +
  coord_flip() +
  theme_minimal() +
  labs(title="Pathway Enrichment Analysis",
       x="Pathway",
       y="Gene Count") +
  theme(
    axis.text.y = element_text(size=9, face="bold", family="serif"),    # Bold y-axis labels
    axis.title.y = element_text(size=12, face="bold", family="serif", colour = "black"),   # Bold y-axis title
    axis.title.x = element_text(size=12, family="serif"),                # X-axis title in Times New Roman
    text = element_text(family="serif")                                  # Everything in Times New Roman
  )
plot(OD1)
library(svglite)
svglite("GSEA_bar_OD1.svg", height = 8, width =8)
print(OD1)
dev.off()

##############
#OD2
library(readxl)
library(ggplot2)
library(dplyr)
setwd("E:/Ssangita_backup/A_Ghosh/Transcriptomics/Deseq2/25032026/")
pathway_data <- read_xlsx("E:/Ssangita_backup/A_Ghosh/2026/KEGG_enrichment_list.xlsx", sheet ="Sorted_OD2")

# Check the data
head(pathway_data)


# Convert Description to factor & order by Count
data <- pathway_data %>% arrange(desc(Count))

# Plot the data with both p-value (color) and FC (bar border)
ggplot(data, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_bar(stat = "identity", color = "black", aes(size = abs(FC))) +  # Bar outline thickness based on FC
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red", name = "P-value") +
  scale_size_continuous(name = "Fold Change", range = c(0.2, 2)) +  # Controls thickness
  labs(title = "Pathway Enrichment Analysis",
       x = "Pathway",
       y = "Gene Count") +
  theme_minimal(base_size = 14)

library(ggplot2)
library(dplyr)

# Assuming your data is loaded in a data frame called df
# df <- read.csv("your_file.csv")

# Sort by Gene Count for better visual representation
df <- data %>%
  arrange(desc(Count))
library(dplyr)
sorted_FC <- data %>% arrange(desc(FC))
# Filter the data frame to include only rows with 'FC' values >= 1.5
final_KEGG_result_sorted_FC <- sorted_FC %>% filter(FC >= 1.5)
# Plotting
OD2 <- ggplot(final_KEGG_result_sorted_FC, aes(x=reorder(Description, Count), y=Count, fill=pvalue)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low="magenta4", high="coral", name="P-value") +
  coord_flip() +
  theme_minimal() +
  labs(title="Pathway Enrichment Analysis",
       x="Pathway",
       y="Gene Count") +
  theme(
    axis.text.y = element_text(size=9, face="bold", family="serif"),    # Bold y-axis labels
    axis.title.y = element_text(size=12, face="bold", family="serif", colour = "black"),   # Bold y-axis title
    axis.title.x = element_text(size=12, family="serif"),                # X-axis title in Times New Roman
    text = element_text(family="serif")                                  # Everything in Times New Roman
  )
plot(OD2)
library(svglite)
svglite("GSEA_bar_OD2.svg", height = 8, width =8)
print(OD2)
dev.off()

