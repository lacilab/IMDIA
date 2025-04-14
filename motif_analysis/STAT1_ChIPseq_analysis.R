###### STAT1 ChIP-Seq analysis

# check whether packages are locally available, if not, install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
  
# load required R packages
library(GenomicRanges)
library(biomaRt)

# inputfiles, change paths here if you wish to rerun
upregulatedgenesfile <- "./STAT1_IFN_upregulated.txt"
statpeaksfile <- "./ENCFF002CTG.bed.gz"

# read in the upregulated genes table and reformat to list
upregulated_genes <- read.delim(file = upregulatedgenesfile, as.is = T)
upregulated_genes_ls <- split(x = t(upregulated_genes), f = 1:ncol(upregulated_genes))
upregulated_genes_ls <- lapply(X = upregulated_genes, FUN = function(x) x[!x == ""])
names(upregulated_genes_ls) <- colnames(upregulated_genes)

# make sure gene sets tested (e.g. WT vs Mut1) are non-overlapping; i.e. remove all WT upregulated genes from Mut upregulated ones
# upregulated_genes_ls_minuswt <- lapply(X = upregulated_genes_ls[-1], FUN = function(x, wtset) setdiff(x = x, y = wtset), wtset = upregulated_genes_ls[[1]])
# upregulated_genes_ls_minuswt[["WT"]] <- upregulated_genes_ls[[1]]
upregulated_genes_ls_minuswt <- upregulated_genes_ls

# load the ENCODE STAT1 ChIP-seq peaks after 30min of IFNg
statpeaks <- read.delim(file = statpeaksfile, header = F, as.is = T)
statpeaks <- GRanges(seqnames = gsub(pattern = "chr", replacement = "", x = statpeaks$V1),
                     ranges = IRanges(start = statpeaks$V2, end = statpeaks$V3))

# get gene loci  (chrom, start, end, strand) from Ensembl biomart
# note that only the first TSS annotated for a gene will be used (i.e. the gene "start" or "end" when on reverse strand)
ensembl <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name"),
               filters = c("external_gene_name", "transcript_gencode_basic", "chromosome_name"), mart = ensembl,
               values = list(external_gene_name = unique(unlist(upregulated_genes_ls)), transcript_gencode_basic = T,
                             chromosome_name = c(1:22, "X", "Y")))

# turn dataframe into GRanges object to make finding promoters and overlaps easy
annotgr <- GRanges(seqnames = annot$chromosome_name, ranges = IRanges(start = annot$start_position, end = annot$end_position),
                   strand = annot$strand, geneid = annot[, -c(1:4)])

# get the promoters for all genes and find overlaps with ChIP-Seq peaks (any overlap will do)
# turn into a named vector where names are genes and elements are TRUE / FALSE if there is an overlapping peak
upreg_promoters <- promoters(x = annotgr, upstream = 2000, downstream = 500)
haspeakv <- setNames(object = overlapsAny(query = upreg_promoters, subject = statpeaks), nm = upreg_promoters$geneid)

# make a small table with peak / no peak counts for all conditions
peakcountstable <- do.call(rbind, lapply(X = upregulated_genes_ls_minuswt, FUN = function(x, haspeakvect) table(haspeakvect[x]), haspeakvect = haspeakv))
print(peakcountstable)

# test 2x2 contingency tables for WT vs each Mutant (i.e. column/row) independence using chi square test
testoutputs <- lapply(X = c("R274W", "R321S", "E411A", "T419R", "N574I"), FUN = function(x, cnttab) chisq.test(x = cnttab[c("WT", x), ]), cnttab = peakcountstable)
names(testoutputs) <- c("R274W", "R321S", "E411A", "T419R", "N574I")

# print outputs
print(testoutputs)
