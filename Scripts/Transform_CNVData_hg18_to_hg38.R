#-------------------------------------------------------------------------------
# The NRC Dataset (GSE85047) was originally mapped to the hg18 reference genome 
# and requires a liftOver to hg19 and then to hg38 for compatibility with 
# updated genomic coordinates and following scripts. 
#-------------------------------------------------------------------------------

# Load libraries
# BiocManager::install("rtracklayer") - V1.62.0
library(rtracklayer)
# BiocManager::install("biomaRt") - V2.58.2
library(biomaRt)
# BiocManager::install("GenomicRanges") - V1.54.1
library(GenomicRanges)
# install.packages("dplyr") - V1.1.4
library(dplyr)
# install.packages("readr") - V2.1.5 
library(readr)

# Set datadir to NRC Dataset
setwd("~/Docker/FilesToPublish/NRC_Dataset/")

# If the file isn't already downloaded, download it manually:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
chain_file <- "hg18ToHg19.over.chain"
chain <- import.chain(chain_file)

# Load Copy Number Variation data 
cnv_data <- read.table("20111106_NRCdata_wavecorrected_en_Dublin_CBS.txt", header = T)
cnv_data$chromosome <- paste0("chr", cnv_data$chromosome)

# Convert to GRanges object
gr_hg18 <- GRanges(
  seqnames = cnv_data$chromosome,
  ranges = IRanges(start = cnv_data$chr_start, end = cnv_data$chr_end), 
  mcols = cnv_data[, c(1,5,6)]) # To keep all the other metadata as well

# Liftover hg18 -> hg19
gr_hg19 <- liftOver(gr_hg18, chain)
# Discarding unchained sequences: chr23, chr24
gr_hg19 <- unlist(gr_hg19)

# If the file isn't already downloaded, download it manually:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg19ToHg38.over.chain.gz
chain_file_2 <- "hg19ToHg38.over.chain"
chain_2 <- import.chain(chain_file_2)

# Liftover hg19 -> hg38
gr_hg38 <- liftOver(gr_hg19, chain_2)
gr_hg38 <- unlist(gr_hg38)

# Into a dataframe
df_hg38 <- data.frame(
  chr = seqnames(gr_hg38),
  start = start(gr_hg38),
  end = end(gr_hg38), 
  mcols(gr_hg38)
)

# Mapping on hg38!
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene annotation
genes_hg38 <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                    filters = "chromosome_name",
                    values = as.character(1:22),  # Only autosomes; add "X", "Y" if needed
                    mart = ensembl)

# Rename columns
colnames(genes_hg38) <- c("chromosome", "start", "end", "gene")

# Ensuring chromosome names are prefixed with "chr".
genes_hg38$chromosome <- paste0("chr", genes_hg38$chromosome)

# Convert hg38 genes into GRanges format
genes_gr_hg38 <- GRanges(
  seqnames = as.character(genes_hg38$chromosome),
  ranges = IRanges(start = genes_hg38$start, end = genes_hg38$end),
  gene = genes_hg38$gene
)
# Convert CNV data into a GRanges object using hg38 annotation
cnv_gr <- GRanges(
  seqnames = as.character(df_hg38$chr),
  ranges = IRanges(start = df_hg38$start, end = df_hg38$end),
  sample = df_hg38$mcols.exp_id,
  ratio = df_hg38$mcols.ratio
)

# Find overlaps
overlaps <- findOverlaps(genes_gr_hg38, cnv_gr)

# Extract matching genes and CNVs
gene_names <- genes_gr_hg38$gene[queryHits(overlaps)]
sample_names <- cnv_gr$sample[subjectHits(overlaps)]
cnv_ratios <- cnv_gr$ratio[subjectHits(overlaps)]

# Create a data frame
cnv_gene_df <- data.frame(
  sample = sample_names,
  gene = gene_names,
  ratio = cnv_ratios
)

# Define CNV thresholds - Based on Depuydt et al. 2018: 10.1038/sdata.2018.240
gain_threshold <- 0.15
loss_threshold <- -0.25

# Summarize CNV status per gene per sample
gene_cnv_status <- cnv_gene_df %>%
  group_by(sample, gene) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE)) %>%
  mutate(conugeal = case_when(
    mean_ratio >= gain_threshold ~ "gain",
    mean_ratio <= loss_threshold ~ "loss",
    TRUE ~ "normal"
  ))

# Output file
write_csv(gene_cnv_status, "gene_level_cnvS_SBK_hg18_hg19_hg38.csv")

