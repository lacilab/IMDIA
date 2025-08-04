# --------------------------- LOAD LIBRARIES --------------------------- #
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sva)
library(pheatmap)
library(clusterProfiler)
library(EnhancedVolcano)

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
sampleInfo$replicate <- factor(sampleInfo$replicate)

#
# --------------------------- DESEQ2 MODEL --------------------------- #
dds <- DESeqDataSetFromMatrix(countData = counts_CB,
                              colData = sampleInfo,
                              design = ~ replicate + mutation + condition)
dds <- DESeq(dds)

# --------------------------- VST for VISUALIZATION --------------------------- #
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# --------------------------- COMBAT for PCA (NOT for DE) --------------------------- #
batch <- sampleInfo$replicate
mod <- model.matrix(~ mutation + condition, data = sampleInfo)
vsd_combat <- ComBat(dat = vsd_mat, batch = batch, mod = mod)

# --------------------------- PCA --------------------------- #
pca <- prcomp(t(vsd_combat))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- sampleInfo$condition
pca_df$mutation <- sampleInfo$mutation
all_shapes <- c(0:25, 32:127)  # 26 default + printable ASCII

ggplot(pca_df, aes(PC1, PC2, color = condition, shape = mutation)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Unstim = "#D3D5D7", IFNa = "#ef3633", IFNg = "#355ba6")) +
  theme_classic() +
  scale_shape_manual(values = all_shapes) +
  labs(title = "PCA Plot of Samples")

# --------------------------- DESEQ2 GROUP MODEL --------------------------- #
sampleInfo$group <- factor(paste0(sampleInfo$mutation, "_", sampleInfo$condition),
                           levels = c("NT_IFNg", setdiff(unique(paste0(sampleInfo$mutation, "_", sampleInfo$condition)), "NT_IFNg")))

dds_grouped <- DESeqDataSetFromMatrix(countData = counts_CB,
                                      colData = sampleInfo,
                                      design = ~ replicate + mutation + condition)
dds_grouped <- DESeq(dds_grouped)

# --------------------------- GENE ANNOTATION --------------------------- #
gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(counts_CB),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)
gene_map_df <- tibble(gene_id = rownames(counts_CB)) %>%
  left_join(gene_annotations, by = c("gene_id" = "ENSEMBL")) %>%
  rename(gene_name = SYMBOL)

# --------------------------- EXPORT DEGs --------------------------- #
ref_group <- "NT_IFNg"
all_groups <- setdiff(levels(sampleInfo$group), ref_group)
dir.create("DEGsNT", showWarnings = FALSE)

for (grp in all_groups) {
  res <- results(dds_grouped, contrast = c("group", grp, ref_group))
  res <- lfcShrink(dds_grouped, contrast = c("group", grp, ref_group), res = res, type = "ashr")

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_annot <- left_join(res_df, gene_map_df, by = "gene_id")

  res_up <- res_annot %>% filter(padj < 0.001, log2FoldChange > 1)
  res_down <- res_annot %>% filter(padj < 0.001, log2FoldChange < -1)

  tag <- gsub("/", "_", grp)
  write.csv(res_annot, file = file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_full.csv")), row.names = FALSE)
  write.csv(res_up, file = file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_up.csv")), row.names = FALSE)
  write.csv(res_down, file = file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_down.csv")), row.names = FALSE)
}

# --------------------------- VOLCANO PLOTS --------------------------- #
dir.create("volcanoNT", showWarnings = FALSE)

for (grp in all_groups) {
  tag <- gsub("/", "_", grp)
  res_file <- file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_full.csv"))
  res <- read.csv(res_file)
  res$gene <- ifelse(is.na(res$gene_name), res$gene_id, res$gene_name)

  p <- EnhancedVolcano(
    res,
    lab = res$gene,
    x = "log2FoldChange",
    y = "padj",
    title = paste0(tag, " vs ", ref_group),
    pCutoff = 0.001,
    FCcutoff = 1
  )
  ggsave(filename = file.path("volcanoNT", paste0("Volcano_", tag, ".pdf")), plot = p, width = 8, height = 6)
}

# --------------------------- FUNCTIONAL ENRICHMENT --------------------------- #
dir.create("enrichmentNT", showWarnings = FALSE)
library(clusterProfiler)
run_enrichment <- function(genes, prefix) {
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  ekegg <- enrichKEGG(gene = genes, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05)

  dotplot(ego, showCategory = 20) + ggtitle(paste0(prefix, " GO:BP"))
  ggsave(paste0("enrichmentNT /", prefix, "_GO_BP.png"), width = 6, height = 5)

  dotplot(ekegg, showCategory = 20) + ggtitle(paste0(prefix, " KEGG"))
  ggsave(paste0("enrichmentNT/", prefix, "_KEGG.png"), width = 6, height = 5)
}

for (grp in all_groups) {
  tag <- gsub("/", "_", grp)
  up_ids <- read.csv(file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_up.csv")))$gene_id
  down_ids <- read.csv(file.path("DEGsNT", paste0(tag, "_vs_", ref_group, "_down.csv")))$gene_id

  run_enrichment(up_ids, paste0(tag, "_up"))
  run_enrichment(down_ids, paste0(tag, "_down"))
}

# --------------------------- HEATMAP OF TOP 50 --------------------------- #
dir.create("plotsNT", showWarnings = FALSE)
example_grp <- all_groups[1]
res_ex <- read.csv(file.path("DEGsNT", paste0(example_grp, "_vs_", ref_group, "_full.csv")))
top50 <- res_ex %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(50) %>% pull(gene_id)

mat_top50 <- assay(vsd)[top50, ]
annotation_col <- sampleInfo[, c("mutation", "condition")]
pheatmap(mat_top50, scale = "row", clustering_method = "ward.D2",
         annotation_col = annotation_col, show_rownames = FALSE,
         filename = "plotsNT/top50_DEG_heatmap.pdf", width = 16, height = 8)

### Z-core Heatmap of Top 50 DEGs
library(pheatmap)
library(dplyr)
library(tibble)

# Ensure sampleInfo$group is defined
sampleInfo$group <- paste0(sampleInfo$mutation, "_", sampleInfo$condition)

# Select one DE comparison
example_grp <- all_groups[1]
res_ex <- read.csv(file.path("DEGs", paste0(example_grp, "_vs_", ref_group, "_full.csv")))

# Select top 50 DEGs
top50 <- res_ex %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(100) %>%
  pull(gene_id)

# Extract VST matrix for top genes
expr_mat <- assay(vsd)[top50, , drop = FALSE]

# Average expression across replicates for each group
grouped_expr <- as.data.frame(expr_mat) %>%
  t() %>%
  as.data.frame() %>%
  mutate(group = sampleInfo$group) %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("group") %>%
  t()  # genes × groups

# Map gene_id to gene_name (with fallback to gene_id if missing)
grouped_expr_df <- as.data.frame(grouped_expr) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_map_df, by = "gene_id") %>%
  mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)) %>%
  select(-gene_id, -gene_name) %>%
  column_to_rownames("label")  # set rownames to gene names (or ID if missing)

# Z-score across groups per gene
z_mat <- t(scale(t(as.matrix(grouped_expr_df))))
z_mat[is.na(z_mat)] <- 0

# Generate annotation for columns
annotation_col <- data.frame(
  mutation = sapply(colnames(z_mat), function(x) strsplit(x, "_")[[1]][1]),
  condition = sapply(colnames(z_mat), function(x) strsplit(x, "_")[[1]][2])
)
rownames(annotation_col) <- colnames(z_mat)

# Save heatmap
pdf("plots/top50_DEG_heatmap_Zscore_grouped_geneNames.pdf", width = 14, height = 14)
pheatmap(z_mat,
         annotation_col = annotation_col,
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         main = paste0("Top 50 DEGs Z-score (gene names) — ", example_grp, " vs ", ref_group))
dev.off()


