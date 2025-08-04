# --------------------------- LOAD LIBRARIES --------------------------- #
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tibble)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# --------------------------- STEP 1: LOAD AND PREPARE DATA --------------------------- #
setwd("/home/luna.kuleuven.be/u0147150/Desktop/PhD/WP2_Overexpression/6_RNAseq/DBD_Muts/Counts/")

# Load counts and metadata
gene_count <- read.delim("gene_count.xls", header = TRUE)
sampleInfo <- read.csv("sampleInfo.csv", header = TRUE, stringsAsFactors = FALSE)

# Assign rownames and slice count matrix
counts <- as.data.frame(gene_count)
rownames(counts) <- counts$gene_id
counts_CB <- counts[ , 2:133]  # adjust if needed
rownames(counts_CB) <- counts$gene_id

# Prepare metadata
rownames(sampleInfo) <- sampleInfo$sample
stopifnot(all(colnames(counts_CB) == rownames(sampleInfo)))

# Factorize variables
sampleInfo$mutation <- factor(sampleInfo$mutation)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c("Unstim", "IFNa", "IFNg"))

# --------------------------- STEP 2: DESEQ2 ANALYSIS --------------------------- #
dds <- DESeqDataSetFromMatrix(
  countData = counts_CB,
  colData = sampleInfo,
  design = ~ mutation + condition
)

# Prefilter genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Relevel condition
dds$condition <- relevel(dds$condition, ref = "Unstim")

# Run DESeq2
dds <- DESeq(dds)

# Get contrasts
res_IFNa <- results(dds, contrast = c("condition", "IFNa", "Unstim"))
res_IFNg <- results(dds, contrast = c("condition", "IFNg", "Unstim"))

# --------------------------- STEP 3: GENE ANNOTATION --------------------------- #
# Retrieve gene symbol mapping from Ensembl IDs
gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(counts_CB),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
) %>%
  dplyr::rename(gene_id = ENSEMBL, gene_name = SYMBOL)

# Create mapping dataframe
gene_map_df <- tibble(gene_id = rownames(counts_CB)) %>%
  left_join(gene_annotations, by = "gene_id")

# --------------------------- STEP 4: PROCESS & EXPORT --------------------------- #
process_and_export <- function(res, gene_map_df, padj_cutoff = 0.05, label = "IFNa") {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # Merge with gene name
  res_df <- left_join(res_df, gene_map_df, by = "gene_id")

  # Disambiguate duplicate symbols
  res_df <- res_df %>%
    group_by(gene_name) %>%
    mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                              paste0(gene_name, "_", row_number()), gene_name)) %>%
    ungroup()

  # Replace rownames with gene_name if available
  rownames(res_df) <- ifelse(!is.na(res_df$gene_name), res_df$gene_name, res_df$gene_id)

  # Save all results
  write.csv(res_df, sprintf("./DESeq2_results_with_GeneNames_%s.csv", label), row.names = TRUE)

  # Save filtered results
  res_filtered <- res_df %>% filter(padj < padj_cutoff)
  write.csv(res_filtered, sprintf("./DESeq2_results_%s_padj%s.csv", label, padj_cutoff), row.names = TRUE)

  # Plots
  plotMA(res, ylim = c(-3, 3), main = paste0("MA Plot - ", label))
  return(res_filtered)
}

# --------------------------- STEP 5: RUN EXPORTS --------------------------- #
# IFNa
process_and_export(res_IFNa, gene_map_df, padj_cutoff = 0.01, label = "IFNa")
process_and_export(res_IFNa, gene_map_df, padj_cutoff = 0.001, label = "IFNa")

# IFNg
process_and_export(res_IFNg, gene_map_df, padj_cutoff = 0.01, label = "IFNg")
process_and_export(res_IFNg, gene_map_df, padj_cutoff = 0.001, label = "IFNg")
























# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(org.Mm.eg.db)  # Mouse gene annotation
library(AnnotationDbi) # Needed for gene ID conversion
library(readxl)
# --------------------------- STEP 1: LOAD AND PREPARE DATA --------------------------- #

setwd("/home/luna.kuleuven.be/u0147150/Desktop/PhD/WP2_Overexpression/6_RNAseq/DBD_Muts/Counts/") 
# Load the gene count data
gene_count <- read.delim("./gene_count.xls", header = TRUE)
colnames(gene_count)
# Load the sample information
sampleInfo <- read.csv("./sampleInfo.csv", header = TRUE, stringsAsFactors = FALSE)
# Check the structure of the data
str(gene_count)
str(sampleInfo)
# Check the first few rows of the data
head(gene_count)
head(sampleInfo)
# Check the column names
colnames(gene_count)
colnames(sampleInfo)
# Check the dimensions of the data
dim(gene_count) 
dim(sampleInfo)
# Check the unique values in the treatment column
unique(sampleInfo$condition)
# Check the unique values in the region column
unique(sampleInfo$mutation)
# Check the unique values in the sampleNames column
# Load count data and sample information
counts <- as.data.frame(gene_count)# Ensure counts is a dataframe

rownames(counts) <- counts$gene_id
counts_CB <- counts[2:133]  # Adjust column indices based on your dataset
rownames(counts_CB) <- counts$gene_id
# Load Sample Info
sampleInfo_CB <- sampleInfo 
rownames(sampleInfo_CB) <- sampleInfo_CB$sample

# Ensure column names match between count matrix and metadata
stopifnot(all(colnames(counts_CB) %in% rownames(sampleInfo_CB)))
stopifnot(all(colnames(counts_CB) == rownames(sampleInfo_CB)))

# --------------------------- STEP 2: DESeq2 ANALYSIS --------------------------- #
sampleInfo$mutation <- factor(sampleInfo$mutation)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c("Unstim", "IFNa", "IFNg"))

# Create DESeq2 object
# Design formula: stimulation condition (IFNa, IFNg, Unstim) and mutant
dds <- DESeqDataSetFromMatrix(
  countData = counts_CB,
  colData = sampleInfo,
  design = ~ mutation + condition
)


# Prefilter low-expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set WT as reference level
dds$condition <- relevel(as.factor(dds$condition), ref = "Unstim")

# Run DESeq2 analysis
dds <- DESeq(dds)
dds <- DESeq(dds)

res_IFNa <- results(dds, contrast = c("condition", "IFNa", "Unstim"))
res_IFNg <- results(dds, contrast = c("condition", "IFNg", "Unstim"))



########### ALLLL ##########
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library(readr)

process_and_export <- function(res, gene_map_df, padj_cutoff = 0.05, label = "IFNa") {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # Merge with gene name
  res_df <- left_join(res_df, gene_map_df, by = "gene_id")

  # Disambiguate duplicated symbols
  res_df <- res_df %>%
    group_by(gene_name) %>%
    mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                              paste0(gene_name, "_", row_number()), gene_name)) %>%
    ungroup()

  # Replace rownames with gene_name if available
  rownames(res_df) <- ifelse(!is.na(res_df$gene_name), res_df$gene_name, res_df$gene_id)

  # Save all results
  write.csv(res_df, sprintf("./DESeq2_results_with_GeneNames_%s.csv", label), row.names = TRUE)

  # Save filtered results
  res_filtered <- res_df %>% filter(padj < padj_cutoff)
  write.csv(res_filtered, sprintf("./DESeq2_results_%s_padj%s.csv", label, padj_cutoff), row.names = TRUE)

  # Plots
  plotMA(res, ylim = c(-3, 3), main = paste0("MA Plot - ", label))
  hist(res$pvalue, breaks = 20, col = "grey", main = paste0("P-value distribution - ", label))
  hist(res$padj, breaks = 20, col = "grey", main = paste0("Adj P-value distribution - ", label))

  return(res_filtered)
}

# Use counts_CB rownames to map gene_id
gene_map_df <- tibble(gene_id = rownames(counts_CB)) %>%
  left_join(gene_annotations, by = "gene_id")  # you must provide gene_annotations

# IFNa vs Unstim
process_and_export(res_ifna, gene_map_df, padj_cutoff = 0.01, label = "IFNa")
process_and_export(res_ifna, gene_map_df, padj_cutoff = 0.001, label = "IFNa")

# IFNg vs Unstim
process_and_export(res_ifng, gene_map_df, padj_cutoff = 0.01, label = "IFNg")
process_and_export(res_ifng, gene_map_df, padj_cutoff = 0.001, label = "IFNg")










# Convert DESeq2 results (res) into a dataframe
deg_data <- as.data.frame(res)

head(deg_data)
# Store rownames (Ensembl IDs) in a new column
deg_data$gene_id <- rownames(deg_data)

# Load your mapping dataframe (Ensure it contains Ensembl IDs and Gene Names)

gene_map_df <- counts %>% dplyr::select(gene_id, gene_name)
# Example: gene_map_df should have columns "ENBL" (Ensembl ID) and "gene_name" (Gene Symbol)

# Merge DESeq2 results with the gene mapping dataframe
deg_data <- left_join(deg_data, gene_map_df, by = "gene_id")
head(deg_data)
deg_data <- deg_data %>%
  group_by(gene_name) %>%
  mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                            paste0(gene_name, "_", row_number()), gene_name)) %>%
  ungroup()
# Replace rownames with gene symbols where available
deg_data <- as.data.frame(deg_data)
rownames(deg_data) <- ifelse(!is.na(deg_data$gene_name), deg_data$gene_name, deg_data$gene_id)

head(deg_data)

# Save the updated DESeq2 results
write.csv(deg_data, file = "./DESeq2_results_with_GeneNames_Condition.csv")

########### padj <0.01 ##############
# Summary with a different padj
res.01 <- results(dds, alpha = 0.01)
summary(res.01)
cook_values <- assays(dds)[["cooks"]]
threshold = 3
high_cooks_per_sample <- colSums(cook_values > threshold, na.rm = TRUE)
barplot(high_cooks_per_sample,
        las = 2,
        main = "Number of Outlier Genes per Sample",
        ylab = paste0("Genes with Cook's D > ", threshold),
        col = "steelblue")


# Convert DESeq2 results (res) into a dataframe
deg_data <- as.data.frame(res.01)

head(deg_data)
# Store rownames (Ensembl IDs) in a new column
deg_data$gene_id <- rownames(deg_data)

# Load your mapping dataframe (Ensure it contains Ensembl IDs and Gene Names)

gene_map_df <- counts %>% dplyr::select(gene_id, gene_name)
# Example: gene_map_df should have columns "ENBL" (Ensembl ID) and "gene_name" (Gene Symbol)

# Merge DESeq2 results with the gene mapping dataframe
deg_data <- left_join(deg_data, gene_map_df, by = "gene_id")
head(deg_data)
deg_data <- deg_data %>%
  group_by(gene_name) %>%
  mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                            paste0(gene_name, "_", row_number()), gene_name)) %>%
  ungroup()
# Replace rownames with gene symbols where available
deg_data <- as.data.frame(deg_data)
rownames(deg_data) <- ifelse(!is.na(deg_data$gene_name), deg_data$gene_name, deg_data$gene_id)

head(deg_data)


# Save the updated DESeq2 results
write.csv(deg_data, file = "./DESeq2_results_with_GeneNames_Condition_p0.01.csv")



####################### padj < 0.001 ###################3


# Summary with a different padj
res.001 <- results(dds, alpha = 0.001)
summary(res.001)
metadata(res.01)$cooksCutoff

# Convert DESeq2 results (res) into a dataframe
deg_data <- as.data.frame(res.001)

head(deg_data)
# Store rownames (Ensembl IDs) in a new column
deg_data$gene_id <- rownames(deg_data)

# Load your mapping dataframe (Ensure it contains Ensembl IDs and Gene Names)

gene_map_df <- counts %>% dplyr::select(gene_id, gene_name)
# Example: gene_map_df should have columns "ENBL" (Ensembl ID) and "gene_name" (Gene Symbol)

# Merge DESeq2 results with the gene mapping dataframe
deg_data <- left_join(deg_data, gene_map_df, by = "gene_id")
head(deg_data)
deg_data <- deg_data %>%
  group_by(gene_name) %>%
  mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                            paste0(gene_name, "_", row_number()), gene_name)) %>%
  ungroup()
# Replace rownames with gene symbols where available
deg_data <- as.data.frame(deg_data)
rownames(deg_data) <- ifelse(!is.na(deg_data$gene_name), deg_data$gene_name, deg_data$gene_id)

head(deg_data)
# Remove the temporary columns used for mapping
deg_data <- deg_data[, !(colnames(deg_data) %in% c("gene_name"))]

# Save the updated DESeq2 results
write.csv(deg_data, file = "./DESeq2_results_with_GeneNames_Condition_p0.001.csv")


# Contrasts: Just in case there are like two treatment conditions instead of one e.g. different time points. Not the case here but will be later (NT, NT-GFP, SAMPLES)
plotMA(res.01, ylim = c(-3,3))
hist(res.01$pvalue, breaks = 20, col = "grey")
hist(res.01$padj, breaks = 20, col = "grey")

plotMA(res.001, ylim = c(-3,3))
hist(res.001$pvalue, breaks = 20, col = "grey")
hist(res.001$padj, breaks = 20, col = "grey")



# --------------------------- STEP 4: VOLCANO PLOT --------------------------- #


### Import the file without cutoffs
DESeq2_results_with_GeneNames_Cerebelum <- read.csv("./DESeq2_results_with_GeneNames_Condition.csv", row.names = 1)
deg_data <- DESeq2_results_with_GeneNames_Cerebelum
# Define thresholds
fc_threshold <- 1
pval_threshold <- 0.001

deg_data$gene_name <- rownames(deg_data)
# Assign colors based on significance
deg_data$color <- "Non-significant"
deg_data$color[deg_data$padj < pval_threshold & deg_data$log2FoldChange > fc_threshold] <- "Upregulated"
deg_data$color[deg_data$padj < pval_threshold & deg_data$log2FoldChange < -fc_threshold] <- "Downregulated"

# Select top 20 upregulated and downregulated DEGs within thresholds
top_upregulated <- deg_data %>%
  filter(padj < pval_threshold & log2FoldChange > fc_threshold) %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

top_downregulated <- deg_data %>%
  filter(padj < pval_threshold & log2FoldChange < -fc_threshold) %>%
  arrange(log2FoldChange) %>%
  head(20)

top_degs <- bind_rows(top_upregulated, top_downregulated)

# Define color palette
color_palette <- c("Non-significant" = "gray70", "Upregulated" = "#ef3633", "Downregulated" = "#355ba6")

# Volcano plot with improved labeling and aesthetics
ggplot(deg_data, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(name = "Differential Expression", values = color_palette) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(data = top_degs, aes(label = gene_name),
                           size = 4, box.padding = 0.5, point.padding = 0.5,
                           max.overlaps = 15, force = 3, segment.color = "gray50") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       subtitle = "Cerebellum: KO vs WT",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 14))


#######################################################################
############################# Heatmap Plot #############################
#######################################################################

vsd <- vst(dds, blind = TRUE)  # or rld <- rlog(dds, blind = TRUE)
vsd_mat <- assay(vsd)
vsd_df <- as.data.frame(vsd_mat)

# Store rownames (Ensembl IDs) in a new column
vsd_df$gene_id <- rownames(vsd_df)

# Load your mapping dataframe (Ensure it contains Ensembl IDs and Gene Names)

gene_map_df <- counts %>% dplyr::select(gene_id, gene_name)
# Example: gene_map_df should have columns "ENBL" (Ensembl ID) and "gene_name" (Gene Symbol)

# Merge DESeq2 results with the gene mapping dataframe
vsd_df <- left_join(vsd_df, gene_map_df, by = "gene_id")
head(vsd_df)
vsd_df <- vsd_df %>%
  group_by(gene_name) %>%
  mutate(gene_name = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                            paste0(gene_name, "_", row_number()), gene_name)) %>%
  ungroup()
# Replace rownames with gene symbols where available
vsd_df <- as.data.frame(vsd_df)
rownames(vsd_df) <- ifelse(!is.na(vsd_df$gene_name), vsd_df$gene_name, vsd_df$gene_id)

head(deg_data)








#deg_filtered <- vsd_df[rownames(vsd_df) %in% vsd_df$gene_name, all_samples]


# Define thresholds
fc_threshold <- 1
pval_threshold <- 0.001


# Select top 50 differentially expressed genes within thresholds
top_degs <- deg_data %>% 
  filter(padj < pval_threshold) %>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  head(50)

# Select top 50 differentially expressed genes within thresholds
top_degs <- deg_data %>% 
  filter(padj < pval_threshold) %>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  head(50)


# Define WT and KO sample groups
wt_samples <- c("CB_01", "CB_02", "CB_03", "CB_07", "CB_08", "CB_09")
ko_samples <- c("CB_04", "CB_05", "CB_06", "CB_10", "CB_11", "CB_12")

# Ensure sample columns exist
all_samples <- c(wt_samples, ko_samples)
sample_cols <- colnames(vsd_df)[colnames(vsd_df) %in% all_samples]
if (length(sample_cols) == 0) {
  stop("No sample expression columns found. Ensure your dataset contains valid sample names.")
}

vsd_df$gene_name <- rownames(vsd_df)
# Extract expression values
deg_filtered <- vsd_df %>% 
  filter(gene_name %in% top_degs$gene_name) %>% 
  dplyr::select(gene_name, all_of(sample_cols))

deg_filtered <- as.data.frame(deg_filtered)
rownames(deg_filtered) <- deg_filtered$gene_name
deg_filtered <- deg_filtered[, -1]  # Remove gene_name column after setting rownames
deg_filtered <- as.matrix(deg_filtered)  # Convert to matrix for heatmap

# Create annotation for column grouping
annotation_col <- data.frame(
  Condition = factor(c(rep("WT", length(wt_samples)), rep("KO", length(ko_samples)))),
  row.names = c(wt_samples, ko_samples)
)

# Define colors for annotation bar
ann_colors <- list(Condition = c(WT = "#D3D5D7", KO = "#CF5C36"))

# Define custom color scale
color_palette <- colorRampPalette(c("#355ba6", "white", "#ef3633"))(100)

library(pheatmap)
# Generate heatmap without hierarchical clustering and with annotation bar
pheatmap(as.matrix(deg_filtered),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Heatmap of Top 50 Differentially Expressed Genes (WT vs KO)",
         color = color_palette)


library(dplyr)

# Filter DEGs with your thresholds
deg_filtered <- deg_data %>%
  filter(padj <= 0.01 & abs(log2FoldChange) > 0.75)

dim(deg_filtered)
library(org.Mm.eg.db)  # Replace with org.Hs.eg.db for human
library(clusterProfiler)

# Convert gene symbols to Entrez IDs
genes_entrez <- bitr(deg_filtered$gene_id, 
                     fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)



ego <- enrichGO(gene          = rownames(deg_filtered),  # your SYMBOL column
                OrgDb         = org.Mm.eg.db,        # use correct organism
                keyType       = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)

head(ego@result$geneID)
# Ensure SYMBOLs are the names
gene_fc <- deg_filtered$log2FoldChange
names(gene_fc) <- rownames(deg_filtered)  # gene should be SYMBOLs

# Initial plot with log2FoldChange coloring
cnetplot(ego,
              showCategory = 10,
              foldChange = gene_fc,
              circular = FALSE,
              colorEdge = F,
              node_label = "all")



goresult <- as.data.frame(ego)
# Generate GO enrichment plot
barplot(ego, showCategory = 30, title = "GO Enrichment Analysis - Biological Process")

################################################################################
go_results <- goresult
# Sort by adjusted p-value (padj)
go_results <- go_results %>% arrange(p.adjust)

# Print the first few rows
print(go_results)

# Convert enrichResult object to a data frame
go_df <- as.data.frame(go_results)

# Export to CSV
write.csv(go_df, file = "go_results.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)

# Ensure Count column exists and is numeric
go_results <- go_results %>% mutate(Count = as.numeric(Count))

library(ggplot2)
library(dplyr)
library(patchwork)  # For combining plots

# Ensure Count column is numeric
go_results <- go_results %>% mutate(Count = as.numeric(Count))
# Convert GeneRatio from "74/1951" format to numeric
# Ensure GeneRatio is calculated correctly from fraction format
go_results <- go_results %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / 
           as.numeric(sub(".*/", "", GeneRatio)))


# Filter by padj < 0.01 and order by Gene Ratio
# Ensure consistent sorting for both plotting & gene extraction
filtered_go <- go_results %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(GeneRatio), p.adjust)  # Ensure stable sorting

# Find global min and max for x-axis scale
gene_ratio_limits <- range(filtered_go$GeneRatio, na.rm = TRUE)

# Set custom limits for color scale (-log10(padj))
padj_limits <- c(0.25, 5)
# Ensure consistent sorting for both plotting & gene extraction
filtered_go <- go_results %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(GeneRatio), p.adjust)  # Ensure stable sorting

# Function to plot per category (top 20 terms)
plot_go <- function(category_name) {
  top_go <- filtered_go %>%
    filter(ONTOLOGY == category_name) %>%
    slice_head(n = 20) %>%  # Select top 20 GO terms **same as gene extraction**
    arrange(GeneRatio)  # Ensure correct order in the plot
  
  # Explicitly set Description as a factor with ordered levels
  top_go$Description <- factor(top_go$Description, levels = top_go$Description)
  
  ggplot(top_go, aes(x = GeneRatio, 
                     y = Description,  
                     color = -log10(p.adjust), size = Count)) +
    geom_point() +
    scale_color_gradient2(low = "#355ba6", mid="#D3D5D7", high = "#CF5C36", limits = padj_limits, oob = scales::squish) +  
    scale_x_continuous(limits = gene_ratio_limits) +  # Set same x-axis scale
    labs(title = paste("Top 20 GO Terms for", category_name),
         x = "Gene Ratio",
         y = "GO Biological Process",
         color = "-log10(adj. p-value)",
         size = "Gene Count") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),  # Angle X-axis labels
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),  # Major grid
      panel.grid.minor = element_line(color = "grey90", linetype = "dotted"),  # Minor grid
      panel.grid.major.x = element_blank()  # Remove horizontal grid if not needed
    )
}

# Generate plots
bp_plot <- plot_go("BP")
cc_plot <- plot_go("CC")
mf_plot <- plot_go("MF")

# Combine plots
pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/GO_Suppl.pdf", width = 12, height = 14)
bp_plot / cc_plot / mf_plot
dev.off()



# Rank all genes by log2FC
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(forcats)

# --------- GSEA Preparation --------- #
gene_ranks <- deg_data$log2FoldChange
names(gene_ranks) <- rownames(deg_data)
gene_ranks <- sort(na.omit(gene_ranks), decreasing = TRUE)

gsea_go <- gseGO(
  geneList      = gene_ranks,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05
)

# --------- Top 20 Terms by NES --------- #
gsea_top <- gsea_go@result %>%
  filter(p.adjust < 0.001 & !is.na(NES)) %>%
  arrange(NES) %>%
  slice_head(n = 20) %>%
  mutate(
    Description = fct_reorder(Description, NES),
    neg_log10_padj = -log10(p.adjust)
  )

# --------- Custom Dotplot --------- #
gsea <- ggplot(gsea_top, aes(x = NES, y = Description)) +
  geom_point(
    aes(size = setSize, fill = neg_log10_padj),
    shape = 21, color = "black", stroke = 0.5
  ) +
  scale_fill_gradientn(
    colors = c("#355ba6", "white", "#ef3633"),
    values = scales::rescale(c(2, 3.5, 5)),
    limits = c(2, 5),
    name = expression(-log[10]~padj)
  ) +
  scale_size(range = c(3, 8), name = "Gene Set Size") +
  labs(
    title = "Top GSEA GO BP Terms (Ranked by NES)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_classic(base_size = 13)
pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/GSEA_GO_BP_Top_20.pdf", width = 8, height = 8)
print(gsea)
dev.off()

# Choose the index (e.g., 1 for top-ranked term)
gs_index <- 1

# Get the gene set ID and title
gs_id <- gsea_top$ID[gs_index]
gs_title <- paste0("NES: ", round(gsea_top$NES[gs_index], 2), " — ", gsea_top$Description[gs_index])

# Plot
p1 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[1], title = gsea_top$Description[1])
p2 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[2], title = gsea_top$Description[2])
p3 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[3], title = gsea_top$Description[3])
p4 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[4], title = gsea_top$Description[4])
p5 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[5], title = gsea_top$Description[5])
p6 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[6], title = gsea_top$Description[6])
p7 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[7], title = gsea_top$Description[7])
p8 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[8], title = gsea_top$Description[8])
p9 <- gseaplot2(gsea_go, geneSetID = gsea_top$ID[9], title = gsea_top$Description[9])





# --------------------- TF STudy ----------------------#
deg_filtered <- deg_data %>%
  filter(padj <= 0.01 & abs(log2FoldChange) > 0.75)
library(enrichR)
dbs <- c("ChEA_2022",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
enriched <- enrichr(rownames(deg_filtered), dbs)

# "ENCODE_Histone_Modifications_2015"

library(tidyr)
library(dplyr)
library(igraph)
library(ggraph)

# Extract top 10 enriched TFs
top_tfs <- enriched[["ChEA_2022"]] %>%
  arrange(desc(Combined.Score)) %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 15)

# Expand TF–target pairs
edge_df <- top_tfs %>%
  select(Term, Genes) %>%
  separate_rows(Genes, sep = ";") %>%
  rename(TF = Term, Target = Genes) %>%
  mutate(
    TF = gsub(" [0-9]+.*", "", TF),
    Target = trimws(Target)
  )

# Create directed graph
tf_net <- graph_from_data_frame(edge_df, directed = TRUE)
tf_list <- unique(edge_df$TF)

# Prepare deg data with uppercase gene symbols
deg_data_upper <- deg_data
rownames(deg_data_upper) <- toupper(rownames(deg_data))

# Assign node attributes
node_df <- data.frame(name = V(tf_net)$name) %>%
  mutate(
    type = ifelse(name %in% tf_list, "TF", "Target"),
    log2FC = ifelse(name %in% rownames(deg_data_upper),
                    deg_data_upper[name, "log2FoldChange"],
                    NA)
  )

V(tf_net)$type <- node_df$type
V(tf_net)$log2FC <- node_df$log2FC
V(tf_net)$degree_out <- degree(tf_net, mode = "out")

# Plot
TF_graph_notOutline_ALL <- ggraph(tf_net, layout = "fr") +
  geom_edge_link(arrow = arrow(length = unit(3, 'mm')), alpha = 0.3) +
  
  # Target nodes
  geom_node_point(
    data = function(x) dplyr::filter(x, type == "Target"),
    aes(size = degree_out, fill = log2FC),
    shape = 21, color = "black", stroke = 0.5
  ) +
  
  # TF nodes
  geom_node_point(
    data = function(x) dplyr::filter(x, type == "TF"),
    aes(size = degree_out + 1),
    shape = 21,
    fill = "#EFC88B",
    color = "black",
    stroke = 0.5
  ) +
  
  # TF labels in bold
  geom_node_text(
    data = function(x) dplyr::filter(x, type == "TF"),
    aes(label = name),
    fontface = "bold",
    size = 3,
    repel = TRUE
  ) +
  
  # Target labels in regular font
  geom_node_text(
    data = function(x) dplyr::filter(x, type == "Target"),
    aes(label = name),
    size = 3,
    repel = TRUE
  ) +
  
  scale_fill_gradient2(
    low = "#355ba6", high = "#ef3633", midpoint = 0,
    limits = c(-2, 2), na.value = "white", name = "log2FC"
  ) +
  scale_size(range = c(3, 8), name = "Outdegree (TF)") +
  theme_void() +
  ggtitle("TF–Target DEG Regulatory Network\nTFs in Beige with Bold Labels")
pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/PublicationReadyTF_graph_notOutline_ALL.pdf", width = 10, height = 8)
print(TF_graph_notOutline_ALL)
dev.off()


# Preprocess for plotting
plot_df <- enriched$`ChEA_2022` %>%
  mutate(
    TF = gsub(" [0-9]+.*", "", Term),  # clean TF names
    logP = -log10(Adjusted.P.value),
    OverlapCount = as.numeric(sub("/.*", "", Overlap)),
    TotalTargets = as.numeric(sub(".*/", "", Overlap)),
    Ratio = OverlapCount / TotalTargets
  ) %>%
  arrange(desc(logP)) %>%
  slice_head(n = 15)  # top 15 enriched TFs

# Dot plot
dotplot_TF_ALL <- ggplot(plot_df, aes(x = logP, y = reorder(TF, logP))) +
  geom_point(aes(size = OverlapCount, color = Combined.Score)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Top Enriched TFs (ChEA_2022)",
    x = "-log10 Adjusted P-value",
    y = "Transcription Factor",
    size = "Overlap Count",
    color = "Combined Score"
  ) +
  theme_classic(base_size = 12)
pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/dotplot_TF_ALL.pdf", width = 10, height = 8)
print(dotplot_TF_ALL)
dev.off()




  
# Heatmap of TF data
  # Required libraries
  library(pheatmap)
  library(dplyr)
  
  # Add gene_name column to vsd_df if not already present
  vsd_df$gene_name <- rownames(vsd_df)
  
  # Convert both to uppercase for matching
  vsd_df$gene_name_upper <- toupper(vsd_df$gene_name)
  node_df$gene_name_upper <- toupper(node_df$name)
  
  # Extract target genes from node_df
  target_genes_upper <- node_df %>%
    filter(type == "Target" & !is.na(log2FC)) %>%
    pull(gene_name_upper)
  
  # Subset vsd_df using uppercased names
  matched_vsd <- vsd_df %>%
    filter(gene_name_upper %in% target_genes_upper)
  
  # Restore original gene names for row labels
  rownames(matched_vsd) <- matched_vsd$gene_name
  
  # Select only expression columns
  matched_vsd <- matched_vsd[, colnames(matched_vsd) %in% c(wt_samples, ko_samples)]
  
  # Ensure numeric matrix
  heatmap_mat <- as.matrix(matched_vsd)
  heatmap_mat <- heatmap_mat[complete.cases(heatmap_mat), ]
  
  # Plot if valid
  if (nrow(heatmap_mat) == 0) stop("Still no matches — check that vsd_df and node_df use the same ID system.")
  # Create annotation for column grouping
  annotation_col <- data.frame(
    Condition = factor(c(rep("WT", length(wt_samples)), rep("KO", length(ko_samples)))),
    row.names = c(wt_samples, ko_samples)
  )
  
  # Define colors for annotation bar
  ann_colors <- list(Condition = c(WT = "#D3D5D7", KO = "#CF5C36"))
  # Define color scale
  color_palette <- colorRampPalette(c("#355ba6", "white", "#ef3633"))(100)
  
  # Generate heatmap
  pheatmap(heatmap_mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = color_palette,
           main = "Unsupervised Clustering of DEG Targets from TF Network")

  deg_data[c("Suz12", "Eed", "Ezh2", "Bmi1", "Rnf2", "Mtf2"), ]
  
  
  
  # -------------- Up and Down Histone modifications ---------------#
  
  up_genes <- rownames(deg_filtered[deg_filtered$log2FoldChange > 0.75 & deg_filtered$padj < 0.01, ])
  down_genes <- rownames(deg_filtered[deg_filtered$log2FoldChange < -0.75 & deg_filtered$padj < 0.01, ])

  library(enrichR)
  dbs <- c("ENCODE_Histone_Modifications_2015")
  
  enrich_up <- enrichr(up_genes, dbs)
  enrich_down <- enrichr(down_genes, dbs)
  
  
  enriched_down_histone <- enrich_down[["ENCODE_Histone_Modifications_2015"]]
  polycomb_genes <- unique(unlist(strsplit(enriched_down_histone$Genes[grep("H3K27me3", enriched_down_histone$Term)], ";")))
  
  # Mark in node_df
  node_df$polycomb_repressed <- ifelse(node_df$name %in% polycomb_genes, TRUE, FALSE)

  
  library(ggplot2)
  library(dplyr)
  
  # Example using enriched histone marks from downregulated genes
  enriched_down_histone <- enrich_down[["ENCODE_Histone_Modifications_2015"]]
  
  # Clean and prepare top 15 entries
  top_histone <- enriched_down_histone %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = 15) %>%
    mutate(Term = factor(Term, levels = rev(Term)))  # reverse for better display
  
  # Plot
  barplot_top_histone_Enrich <- ggplot(top_histone, aes(x = Term, y = -log10(Adjusted.P.value), fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "#6baed6", high = "#2171b5") +
    labs(
      title = "Top Enriched Histone Marks in Downregulated Genes",
      x = "Histone Modification (ENCODE)",
      y = "-log10 Adjusted P-value",
      fill = "Combined Score"
    ) +
    theme_minimal(base_size = 12)
  
  
  
  

  
  
  
  
  library(enrichR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  
  ###--- Assumes you already have:
  # - `enriched_down` from enrichr()
  # - `deg_filtered`: filtered DEGs
  # - `deg_data`: full DEG table with log2FC
  # - `tf_net`: igraph object of TF–target network
  # - `node_df`: data.frame of node attributes

  deg_filtered <- deg_data %>%
  filter(padj <= 0.01 & abs(log2FoldChange) > 0.75)

  ### 1. Extract H3K27me3 CNS-relevant terms
  h3k27_terms <- enrich_down[["ENCODE_Histone_Modifications_2015"]] %>%
    filter(str_detect(Term, "H3K27me3"),
           str_detect(Term, "(?i)astrocyte|SK-N-SH|brain|ES|hESC|neur|mono"))
  
  # Flatten gene list
  polycomb_genes <- h3k27_terms$Genes %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    toupper() %>%
    unique()
  
  ### 2. Flag nodes in node_df
  canonical_polycomb_targets <- toupper(c("HOXA1", "HOXA2", "HOXA3", "HOXD3", "MEIS1"))
  
  node_df <- node_df %>%
    mutate(
      polycomb_flag = name %in% polycomb_genes,
      canonical_polycomb = name %in% canonical_polycomb_targets,
      color_border = ifelse(polycomb_flag, "purple", "black")
    )
  
  # Update graph with new attributes
  V(tf_net)$type <- node_df$type
  V(tf_net)$log2FC <- node_df$log2FC
  V(tf_net)$degree_out <- degree(tf_net, mode = "out")
  V(tf_net)$canonical_polycomb <- node_df$canonical_polycomb
  V(tf_net)$color_border <- node_df$color_border
  
  ### 3. Plot
  ggraph(tf_net, layout = "fr") +
    geom_edge_link(arrow = arrow(length = unit(3, 'mm')), alpha = 0.3) +
    
    # Target nodes
    geom_node_point(
      data = function(x) filter(x, type == "Target"),
      aes(size = degree_out, fill = log2FC, color = I(color_border)),
      shape = 21, stroke = 0.7
    ) +
    
    # TF nodes
    geom_node_point(
      data = function(x) filter(x, type == "TF"),
      aes(size = degree_out + 1),
      shape = 21, fill = "#EFC88B", color = "black", stroke = 0.5
    ) +
    
    # Canonical Polycomb target labels
    geom_node_text(
      data = function(x) filter(x, canonical_polycomb),
      aes(label = name),
      fontface = "bold", size = 3.5, color = "red", repel = TRUE
    ) +
    
    # TF labels
    geom_node_text(
      data = function(x) filter(x, !canonical_polycomb & type == "TF"),
      aes(label = name),
      fontface = "bold", size = 3, repel = TRUE
    ) +
    
    # Target labels
    geom_node_text(
      data = function(x) filter(x, !canonical_polycomb & type == "Target"),
      aes(label = name),
      size = 2.7, repel = TRUE
    ) +
    
    scale_fill_gradient2(
      low = "#355ba6", high = "#ef3633", midpoint = 0,
      limits = c(-2, 2), na.value = "white", name = "log2FC"
    ) +
    scale_size(range = c(3, 8), name = "Outdegree (TF)") +
    theme_void() +
    ggtitle("TF–Target DEG Network with Polycomb Repression Overlay")
  
  

  ## Coloring for synapse GO term
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)

# Define GO term of interest: Synapse
go_term <- "GO:0045202"

# Get all mappings for GO:0045202
synapse_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_term,
  keytype = "GO",
  columns = c("SYMBOL")
) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)

go_synapse_terms <- c(
  "GO:0045202", "GO:0098793", "GO:0098794", "GO:0008021", "GO:0007268",
  "GO:0007269", "GO:0097485", "GO:0048487", "GO:0099536", "GO:0048856",
  "GO:0097110", "GO:0014069", "GO:0014068", "GO:0098978", "GO:0098982"
)

library(org.Mm.eg.db)
library(AnnotationDbi)

# Get gene symbols annotated to any of the synaptic GO terms
expanded_synapse_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_synapse_terms,
  keytype = "GO",
  columns = "SYMBOL"
) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)


node_df <- node_df %>%
  mutate(
    synaptic_flag = toupper(name) %in% toupper(expanded_synapse_genes),
    color_border = ifelse(synaptic_flag, "green", "black")
  )

# Ensure attributes are assigned to graph from node_df
V(tf_net)$type <- node_df$type
V(tf_net)$log2FC <- node_df$log2FC
V(tf_net)$degree_out <- degree(tf_net, mode = "out")
V(tf_net)$synaptic_flag <- node_df$synaptic_flag
V(tf_net)$color_border <- node_df$color_border

  
  ### 3. Plot
library(ggraph)
library(ggplot2)

TF_graph_Synap_ALL <- ggraph(tf_net, layout = "fr") +
  # Edges
  geom_edge_link(arrow = arrow(length = unit(3, 'mm')), alpha = 0.3) +
  
  # Target nodes
  geom_node_point(
    data = function(x) filter(x, type == "Target"),
    aes(size = degree_out, fill = log2FC, color = I(color_border)),
    shape = 21, stroke = 0.7
  ) +
  
  # TF nodes
  geom_node_point(
    data = function(x) filter(x, type == "TF"),
    aes(size = degree_out + 1),
    shape = 21, fill = "#EFC88B", color = "black", stroke = 0.5
  ) +
  
  # Synaptic gene labels
  geom_node_text(
    data = function(x) filter(x, synaptic_flag),
    aes(label = name),
    fontface = "bold", size = 3.5, color = "darkgreen", repel = TRUE
  ) +
  
  # Other TF labels
  geom_node_text(
    data = function(x) filter(x, !synaptic_flag & type == "TF"),
    aes(label = name),
    fontface = "bold", size = 3, repel = TRUE
  ) +
  
  # Other target labels
  geom_node_text(
    data = function(x) filter(x, !synaptic_flag & type == "Target"),
    aes(label = name),
    size = 2.7, repel = TRUE
  ) +
  
  # Scales and theme
  scale_fill_gradient2(
    low = "#355ba6", high = "#ef3633", midpoint = 0,
    limits = c(-2, 2), na.value = "white", name = "log2FC"
  ) +
  scale_size(range = c(3, 8), name = "Outdegree (TF)") +
  theme_void(base_size = 13) +
  ggtitle("TF–Target DEG Network with Synaptic Gene Highlight")

pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/TF_graph_Synap_ALL.pdf", width = 10, height = 8)
print(TF_graph_Synap_ALL)
dev.off()

# --------------------- TF STudy Split by Up/Down ----------------------#
# -------------------- Load Required Packages -------------------- #
library(enrichR)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyr)
library(dplyr)
library(igraph)
library(ggraph)
library(ggnewscale)
library(ggplot2)

# -------------------- DEG Filtering -------------------- #
lfc_thresh <- 0.75
padj_thresh <- 0.001

deg_up <- deg_data %>%
  filter(padj <= padj_thresh & log2FoldChange > lfc_thresh)

deg_down <- deg_data %>%
  filter(padj <= padj_thresh & log2FoldChange < -lfc_thresh)

# -------------------- TF Enrichment -------------------- #
dbs <- c("ChEA_2022", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

enriched_up <- enrichr(rownames(deg_up), dbs)
enriched_down <- enrichr(rownames(deg_down), dbs)

# -------------------- Synaptic Gene Set -------------------- #
go_terms_synapse <- c(
  "GO:0045202",  # Synapse
  "GO:0098793",  # Presynapse
  "GO:0098794",  # Postsynapse
  "GO:0007268",  # Synaptic transmission
  "GO:0050804",  # Modulation of synaptic transmission
  "GO:0007411",   # Axon guidance
  "GO:0051963", # synapse assembly
  "GO:0051965"
)

go_synapse_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_terms_synapse,
  keytype = "GO",
  columns = c("SYMBOL", "GO")
)

expanded_synapse_genes <- toupper(unique(go_synapse_genes$SYMBOL))

# -------------------- TF–Target Network Builder -------------------- #
build_tf_network <- function(enriched_df, deg_subset, top_n = 15) {
  top_tfs <- enriched_df[["ChEA_2022"]] %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = top_n)

  edge_df <- top_tfs %>%
    select(Term, Genes) %>%
    separate_rows(Genes, sep = ";") %>%
    rename(TF = Term, Target = Genes) %>%
    mutate(
      TF = gsub(" [0-9]+.*", "", TF),
      Target = trimws(Target)
    )

  tf_list <- unique(edge_df$TF)
  deg_upper <- deg_subset
  rownames(deg_upper) <- toupper(rownames(deg_upper))

  tf_net <- graph_from_data_frame(edge_df, directed = TRUE)

  node_df <- data.frame(name = V(tf_net)$name) %>%
    mutate(
      name_upper = toupper(name),
      type = ifelse(name %in% tf_list, "TF", "Target"),
      log2FC = ifelse(name_upper %in% rownames(deg_upper),
                      deg_upper[name_upper, "log2FoldChange"],
                      NA),
      synaptic_flag = name_upper %in% expanded_synapse_genes,
      color_border = ifelse(synaptic_flag, "green", "black")
    )

  V(tf_net)$type <- node_df$type
  V(tf_net)$log2FC <- node_df$log2FC
  V(tf_net)$degree_out <- degree(tf_net, mode = "out")
  V(tf_net)$synaptic_flag <- node_df$synaptic_flag
  V(tf_net)$color_border <- node_df$color_border

  return(list(graph = tf_net, nodes = node_df, edges = edge_df, top_tfs = top_tfs))
}

# -------------------- Network Plot Function -------------------- #
plot_tf_network_synaptic <- function(tf_net, title = "", max_outdegree = NULL) {
  # If max_outdegree not provided, compute from this graph
  if (is.null(max_outdegree)) {
    max_outdegree <- max(V(tf_net)$degree_out, na.rm = TRUE)
  }

  ggraph(tf_net, layout = "fr") +
    # Edges
    geom_edge_link(arrow = arrow(length = unit(3, 'mm')), alpha = 0.3) +

    # Target nodes
    geom_node_point(
      data = function(x) filter(x, type == "Target"),
      aes(size = degree_out, fill = log2FC, color = I(color_border)),
      shape = 21, stroke = 0.8
    ) +

    # TF nodes
    new_scale("color") +
    geom_node_point(
      data = function(x) filter(x, type == "TF"),
      aes(size = degree_out + 1),
      shape = 21, fill = "#EFC88B", color = "black", stroke = 0.5
    ) +

    # Synaptic labels
    geom_node_text(
      data = function(x) filter(x, synaptic_flag),
      aes(label = name),
      fontface = "bold", size = 3.5, color = "darkgreen", repel = TRUE
    ) +

    # TF labels
    geom_node_text(
      data = function(x) filter(x, !synaptic_flag & type == "TF"),
      aes(label = name),
      fontface = "bold", size = 3, repel = TRUE
    ) +

    # Other target labels
    geom_node_text(
      data = function(x) filter(x, !synaptic_flag & type == "Target"),
      aes(label = name),
      size = 2.7, repel = TRUE
    ) +

    # Fill and size scales
    scale_fill_gradient2(
      low = "#355ba6", high = "#ef3633", midpoint = 0,
      limits = c(-2, 2), na.value = "white", name = "log2FC"
    ) +
    scale_size(
      limits = c(0, max_outdegree),
      range = c(3, 8),
      name = "Outdegree (TF)"
    ) +

    theme_void(base_size = 13) +
    ggtitle(title)
}

# -------------------- Build & Plot Networks -------------------- #
up_result <- build_tf_network(enriched_up, deg_up, top_n = 10)
down_result <- build_tf_network(enriched_down, deg_down, top_n = 10)



up_graph_synap <- plot_tf_network_synaptic(up_result$graph,
                         title = "Upregulated TF–Target Network with Synaptic Highlight",
                         max_outdegree = 100)
down_graph_synap <- plot_tf_network_synaptic(down_result$graph,
                         title = "Downregulated TF–Target Network with Synaptic Highlight",
                         max_outdegree = 100)
TF_graph_Synap_Split <- up_graph_synap | down_graph_synap

pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/TF_graph_Synap_Split.pdf", width = 20, height = 8)
print(TF_graph_Synap_Split)
dev.off()






non_synaptic_down <- down_result$nodes %>%
  filter(type == "Target", log2FC < -0.75, !synaptic_flag) %>%
  pull(name) %>%
  unique()

non_synaptic_down <- str_to_title(non_synaptic_down)  # From stringr

go_enrich <- enrichGO(
  gene         = non_synaptic_down,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable     = TRUE
)

head(go_enrich@result, 20) 


# --------------------- TF STudy Split by Up/Down ----------------------#
# -------------------- Load Required Packages -------------------- #
library(enrichR)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyr)
library(dplyr)
library(igraph)
library(ggraph)
library(ggnewscale)
library(ggplot2)

# -------------------- DEG Filtering -------------------- #
lfc_thresh <- 0.75
padj_thresh <- 0.001

deg_up <- deg_data %>%
  filter(padj <= padj_thresh & log2FoldChange > lfc_thresh)

deg_down <- deg_data %>%
  filter(padj <= padj_thresh & log2FoldChange < -lfc_thresh)

# -------------------- TF Enrichment -------------------- #
dbs <- c("ChEA_2022", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

enriched_up <- enrichr(rownames(deg_up), dbs)
enriched_down <- enrichr(rownames(deg_down), dbs)

# -------------------- Synaptic Gene Set -------------------- #
go_terms_synapse <- c(
  "GO:0045202",  # Synapse
  "GO:0098793",  # Presynapse
  "GO:0098794",  # Postsynapse
  "GO:0007268",  # Synaptic transmission
  "GO:0050804",  # Modulation of synaptic transmission
  "GO:0007411",   # Axon guidance
  "GO:0051963", # synapse assembly
  "GO:0051965"
)

go_synapse_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = go_terms_synapse,
  keytype = "GO",
  columns = c("SYMBOL", "GO")
)

expanded_synapse_genes <- toupper(unique(go_synapse_genes$SYMBOL))

# -------------------- TF–Target Network Builder -------------------- #
build_tf_network <- function(enriched_df, deg_subset, top_n = 15) {
  top_tfs <- enriched_df[["ChEA_2022"]] %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = top_n)

  edge_df <- top_tfs %>%
    select(Term, Genes) %>%
    separate_rows(Genes, sep = ";") %>%
    rename(TF = Term, Target = Genes) %>%
    mutate(
      TF = gsub(" [0-9]+.*", "", TF),
      Target = trimws(Target)
    )

  tf_list <- unique(edge_df$TF)
  deg_upper <- deg_subset
  rownames(deg_upper) <- toupper(rownames(deg_upper))

  tf_net <- graph_from_data_frame(edge_df, directed = TRUE)

  node_df <- data.frame(name = V(tf_net)$name) %>%
    mutate(
      name_upper = toupper(name),
      type = ifelse(name %in% tf_list, "TF", "Target"),
      log2FC = ifelse(name_upper %in% rownames(deg_upper),
                      deg_upper[name_upper, "log2FoldChange"],
                      NA),
      synaptic_flag = name_upper %in% expanded_synapse_genes,
      color_border = ifelse(synaptic_flag, "green", "black")
    )

  V(tf_net)$type <- node_df$type
  V(tf_net)$log2FC <- node_df$log2FC
  V(tf_net)$degree_out <- degree(tf_net, mode = "out")
  V(tf_net)$synaptic_flag <- node_df$synaptic_flag
  V(tf_net)$color_border <- node_df$color_border

  return(list(graph = tf_net, nodes = node_df, edges = edge_df, top_tfs = top_tfs))
}

# -------------------- Network Plot Function -------------------- #
plot_tf_network_synaptic <- function(tf_net, title = "", max_outdegree = NULL) {
  # If max_outdegree not provided, compute from this graph
  if (is.null(max_outdegree)) {
    max_outdegree <- max(V(tf_net)$degree_out, na.rm = TRUE)
  }

  ggraph(tf_net, layout = "fr") +
    # Edges
    geom_edge_link(arrow = arrow(length = unit(3, 'mm')), alpha = 0.3) +

    # Target nodes
    geom_node_point(
      data = function(x) filter(x, type == "Target"),
      aes(size = degree_out, fill = log2FC, color = I(color_border)),
      shape = 21, stroke = 0.8
    ) +

    # TF nodes
    new_scale("color") +
    geom_node_point(
      data = function(x) filter(x, type == "TF"),
      aes(size = degree_out + 1),
      shape = 21, fill = "#EFC88B", color = "black", stroke = 0.5
    ) +

    # Synaptic labels
    geom_node_text(
      data = function(x) filter(x, synaptic_flag),
      aes(label = name),
      fontface = "bold", size = 3.5, color = "darkgreen", repel = TRUE
    ) +

    # TF labels
    geom_node_text(
      data = function(x) filter(x, !synaptic_flag & type == "TF"),
      aes(label = name),
      fontface = "bold", size = 3, repel = TRUE
    ) +

    # Other target labels
    geom_node_text(
      data = function(x) filter(x, !synaptic_flag & type == "Target"),
      aes(label = name),
      size = 2.7, repel = TRUE
    ) +

    # Fill and size scales
    scale_fill_gradient2(
      low = "#355ba6", high = "#ef3633", midpoint = 0,
      limits = c(-2, 2), na.value = "white", name = "log2FC"
    ) +
    scale_size(
      limits = c(0, max_outdegree),
      range = c(3, 8),
      name = "Outdegree (TF)"
    ) +

    theme_void(base_size = 13) +
    ggtitle(title)
}

# -------------------- Build & Plot Networks -------------------- #
up_result <- build_tf_network(enriched_up, deg_up, top_n = 20)
down_result <- build_tf_network(enriched_down, deg_down, top_n = 20)



up_graph_synap <- plot_tf_network_synaptic(up_result$graph,
                         title = "Upregulated TF–Target Network with Synaptic Highlight",
                         max_outdegree = 75)
down_graph_synap <- plot_tf_network_synaptic(down_result$graph,
                         title = "Downregulated TF–Target Network with Synaptic Highlight",
                         max_outdegree = 75)
TF_graph_Synap_Split <- up_graph_synap | down_graph_synap

pdf("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/TF_graph_Synap_Split.pdf", width = 20, height = 8)
print(TF_graph_Synap_Split)
dev.off()





# --------------------- TF–Target Network via ChEA_2022 ---------------------- #

# --------------------- TF–Target Network via ChEA_2022 ---------------------- #

# Load required libraries
library(enrichR)
library(dplyr)
library(tidyr)
library(tibble)
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)
library(patchwork)

# -------------------- DEG Filtering -------------------- #
lfc_thresh <- 0.75
padj_thresh <- 0.001

deg_data <- deg_data %>%
  as.data.frame() %>%
  rownames_to_column("GENE") %>%
  mutate(GENE = toupper(GENE)) %>%
  column_to_rownames("GENE")

up_genes <- rownames(deg_data[deg_data$padj < padj_thresh & deg_data$log2FoldChange > lfc_thresh, ])
down_genes <- rownames(deg_data[deg_data$padj < padj_thresh & deg_data$log2FoldChange < -lfc_thresh, ])

# -------------------- Enrichment -------------------- #
dbs <- c("ChEA_2022")
enriched_up <- enrichr(up_genes, dbs)
trrust_up <- enriched_up[["ChEA_2022"]] %>% mutate(Direction = "Upregulated")
enriched_down <- enrichr(down_genes, dbs)
trrust_down <- enriched_down[["ChEA_2022"]] %>% mutate(Direction = "Downregulated")

# -------------------- Optional Functional Gene Set (e.g., Synapse) -------------------- #
# expanded_synapse_genes <- c("SYT1", "SNAP25", ...)  # or NULL if not used

# -------------------- Network Plotting Function -------------------- #
plot_trrust_network_corrected <- function(tf_table,
                                          deg_data,
                                          direction_color = "#ef3633",
                                          top_n = 10,
                                          expanded_functional_genes = NULL,
                                          plot_title = NULL) {
  # ---- Prepare DEG data (ensure uppercase gene symbols) ----
  deg_data <- deg_data %>%
    as.data.frame() %>%
    rownames_to_column("GENE") %>%
    mutate(GENE = toupper(GENE)) %>%
    column_to_rownames("GENE")

  # ---- Top TF selection: collapse ChIP variants to TF symbols ----
  tf_table <- tf_table %>%
    mutate(Term = toupper(gsub(" [0-9]+.*", "", Term))) %>%  # Clean TF names
    group_by(Term) %>%
    slice_min(Adjusted.P.value, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = top_n) %>%
    mutate(Genes = toupper(Genes))

  # ---- Build edges ----
  edges <- tf_table %>%
    select(Term, Genes, Direction) %>%
    separate_rows(Genes, sep = ";") %>%
    mutate(
      Term = trimws(Term),
      Gene = trimws(Genes),
      TF = Term
    ) %>%
    select(TF, Gene, Direction) %>%
    mutate(Gene = toupper(Gene), TF = toupper(TF)) %>%
    distinct()

  # ---- Build node metadata ----
  all_nodes <- unique(c(edges$TF, edges$Gene))

  deg_df <- deg_data %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    select(name, log2FoldChange)

  node_df <- tibble(name = all_nodes) %>%
    mutate(
      type = ifelse(name %in% edges$TF, "TF", "Target"),
      synaptic_flag = if (!is.null(expanded_functional_genes)) name %in% toupper(expanded_functional_genes) else FALSE
    ) %>%
    left_join(deg_df, by = c("name"))  # log2FC is NA if not in DEGs

  # ---- Build graph ----
  g <- graph_from_data_frame(edges, vertices = node_df, directed = TRUE)

  # ---- Plot ----
  ggraph(g, layout = "fr") +
    geom_edge_link(color = direction_color, alpha = 0.6,
                   arrow = arrow(length = unit(2, 'mm'))) +
    geom_node_point(
      aes(shape = type, fill = log2FoldChange,
          stroke = ifelse(type == "TF", 1.5, 0.3)),
      size = 5,
      color = ifelse(V(g)$synaptic_flag, "darkgreen", "black")
    ) +
    geom_node_text(
      aes(label = name, fontface = ifelse(type == "TF", "bold", "plain")),
      repel = TRUE, size = 3,
      color = ifelse(V(g)$synaptic_flag, "darkgreen", "black")
    ) +
    scale_shape_manual(values = c(TF = 21, Target = 21)) +
    scale_fill_gradient2(
      low = "#355ba6", mid = "white", high = "#ef3633",
      midpoint = 0, limits = c(-2, 2), oob = squish, na.value = "grey90",
      name = "log2FC"
    ) +
    labs(
      title = if (!is.null(plot_title)) plot_title else paste("ChEA_2022 TF–Target Network:", tf_table$Direction[1]),
      subtitle = paste("Top", top_n, "enriched TFs")
    ) +
    theme_void()
}


# -------------------- Generate and Save Plots -------------------- #

p_up <- plot_trrust_network_corrected(
  tf_table = trrust_up,
  deg_data = deg_data,
  direction_color = "#ef3633",
  top_n = 10,
  expanded_functional_genes = NULL,  # or your expanded_synapse_genes
  plot_title = "Upregulated TF–Target Network"
)

p_down <- plot_trrust_network_corrected(
  tf_table = trrust_down,
  deg_data = deg_data,
  direction_color = "#355ba6",
  top_n = 10,
  expanded_functional_genes = NULL,  # or your expanded_synapse_genes
  plot_title = "Downregulated TF–Target Network"
)

final_plot <- p_up + p_down + plot_layout(ncol = 1)

ggsave("/home/luna.kuleuven.be/u0147150/Desktop/PhD/Other/Laura/Analysis/TF_graph_ChEA2022_Split.pdf",
       plot = final_plot, width = 10, height = 14)



