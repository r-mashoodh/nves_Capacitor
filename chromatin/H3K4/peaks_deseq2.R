setwd("~/Dropbox/Lab/Projects/Beetles/Genome_Reassembly/H3K4me3/peaks/")  # update this path
library(DESeq2)
library(tidyverse)
library(GenomicRanges)

# Read counts
counts <- read.table("larvae.H3K4.counts") %>%
  select(-c(V5, V6, V7, V8, V9, V10)) #remove narrowPeak file cols

colnames(counts) <- c("seqnames",
                      "start",
                      "end",
                      "peak",
                      "CT1_FC", "CT1_NC", "CT2_FC", "CT2_NC",
                      "CT3_FC", "CT3_NC", "CT4_FC", "CT4_NC")
head(counts)
nrow(counts)

# Remove chrMT and unplaced scaffolds
keep_chr <- !grepl("^(chrMT|scaffold_)", counts$seqnames)
counts <- counts[keep_chr, ]
nrow(counts)

# Start DESeq2
countData <- counts %>%
  select(4:12) %>%
  column_to_rownames(var = "peak")

countData <- as.matrix(countData)

metaData <- data.frame(ID = c("CT1_FC", "CT1_NC", "CT2_FC", "CT2_NC",
                              "CT3_FC", "CT3_NC", "CT4_FC", "CT4_NC"))
metaData$Cond <- as.factor(substr(metaData$ID, 5, 6))
metaData$Rep  <- as.factor(c("r1", "r1", "r2", "r2", "r3", "r3", "r4", "r4"))

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ Rep + Cond)

# Remove low count peaks
#dds <- dds[rowMeans(counts(dds)) > 100, ]
nrow(dds)

# Transformation for QC plots
rld <- rlog(dds, blind = FALSE)

plotPCA(rld, intgroup = "Cond", ntop = 10000) +
  geom_text(aes(label = rld$ID), nudge_y = 1, size = 3) +
  ggtitle("H3K4me3 - PCA by condition")

plotPCA(rld, intgroup = "Rep", ntop = nrow(rld)) +
  geom_text(aes(label = rld$ID), nudge_y = 1, size = 3) +
  ggtitle("H3K4me3 - PCA by replicate")

colSums(counts(dds))

# Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = "Cond_NC_vs_FC")

length(which(res$padj < 0.05))
# [1] 1174

options(scipen = 999)
plotMA(res, ylim = c(-4, 4))
title("H3K4me3 Genrich peaks - NC vs FC")

res

write_tsv(data.frame(res), "deseq_da_peaks.tsv")


sig_peaks <- data.frame(res) %>% 
  filter(padj < 0.05)

counts %>% 
  filter(peak %in% rownames(sig_peaks)) %>% 
  select(seqnames,start,end,peak)
  
DE <- read_tsv("../../NewRNAseq/results/Care_Deseq2.tsv")
head(DE)

library(ChIPseeker)
library(GenomicFeatures)
library(tidyverse)
library(GenomicRanges)

# Build TxDb from your Ni2 GTF
txdb <- makeTxDbFromGFF(
  "~/new_hic_genome/Nvi2.liftoff.mitohifi.final.ok.genename.for.cellranger.gtf",
  format = "gtf"
)

# Make GRanges from significant peaks
sig_gr <- counts %>%
  filter(peak %in% rownames(sig_peaks)) %>%
  dplyr::select(seqnames, start, end, peak) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Annotate peaks to nearest gene
anno <- annotatePeak(sig_gr,
                     tssRegion = c(-3000, 3000),
                     TxDb = txdb,
                     annoDb = NULL)

anno_df <- as.data.frame(anno)
head(anno_df)

# Join with DE results
integrated <- anno_df %>%
  select(seqnames, start, end, peak, annotation, distanceToTSS, geneId) %>%
  left_join(DE, by = c("geneId" = "gene_id"))

head(integrated)

# How many significant peaks are near DE genes?
integrated %>%
  filter(padj < 0.05) %>%
  nrow()

# Split by direction
integrated %>%
  filter(padj < 0.05) %>%
  count(log2FoldChange > 0)

# Plot - H3K4me3 fold change vs RNA-seq fold change
integrated %>%
  filter(!is.na(log2FoldChange)) %>%
  ggplot(aes(x = log2FoldChange, y = distanceToTSS)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "RNA-seq log2FC (NC vs FC)",
       y = "Distance to TSS",
       title = "H3K4me3 peaks vs differential expression") +
  theme_bw()

# More useful - just show the DE genes near peaks
integrated %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  select(peak, seqnames, start, end, annotation, distanceToTSS, 
         gene_name, log2FoldChange, padj) %>%
  print(n = 20)

# Write out
write_tsv(integrated, "H3K4_peaks_DE_integrated.tsv")

# Add direction labels
integrated <- integrated %>%
  mutate(
    H3K4_direction = case_when(
      peak %in% rownames(sig_peaks[sig_peaks$log2FoldChange > 0,]) ~ "NC_enriched",
      peak %in% rownames(sig_peaks[sig_peaks$log2FoldChange < 0,]) ~ "FC_enriched",
      TRUE ~ "not_sig"
    ),
    RNA_direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "NC_up",
      padj < 0.05 & log2FoldChange < 0 ~ "FC_up",
      TRUE ~ "not_DE"
    ),
    concordant = case_when(
      H3K4_direction == "NC_enriched" & RNA_direction == "NC_up" ~ "concordant",
      H3K4_direction == "FC_enriched" & RNA_direction == "FC_up" ~ "concordant",
      H3K4_direction == "NC_enriched" & RNA_direction == "FC_up" ~ "discordant",
      H3K4_direction == "FC_enriched" & RNA_direction == "NC_up" ~ "discordant",
      TRUE ~ "not_sig"
    )
  )

# Summary
integrated %>% count(H3K4_direction, RNA_direction, concordant)

# Concordant peaks with gene names
integrated %>%
  filter(concordant == "concordant") %>%
  arrange(padj) %>%
  select(peak, seqnames, start, end, H3K4_direction, 
         gene_name, log2FoldChange, padj)

# Plot
integrated %>%
  filter(!is.na(log2FoldChange), H3K4_direction != "not_sig") %>%
  ggplot(aes(x = log2FoldChange, fill = H3K4_direction)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~H3K4_direction) +
  labs(x = "RNA-seq log2FC (NC vs FC)",
       title = "H3K4me3 direction vs RNA-seq expression") +
  theme_bw()

integrated %>%
  filter(abs(distanceToTSS) < 1500) %>%
  count(H3K4_direction, RNA_direction, concordant)

# NC_enriched peaks: concordant vs discordant (excluding not_DE)
nc_mat <- matrix(c(213, 148, 68, 156), nrow = 2,
                 dimnames = list(c("NC_enriched", "FC_enriched"),
                                 c("NC_up", "FC_up")))
nc_mat
fisher.test(nc_mat)

fc_mat <- matrix(c(156, 68, 148, 213), nrow = 2,
                 dimnames = list(c("FC_enriched", "NC_enriched"),
                                 c("FC_up", "NC_up")))
fc_mat
fisher.test(fc_mat)

integrated %>%
  filter(abs(distanceToTSS) < 1500, !is.na(log2FoldChange)) %>%
  ggplot(aes(x = H3K4_direction, y = log2FoldChange, fill = H3K4_direction)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("NC_enriched" = "#2166ac", "FC_enriched" = "#d73027", "not_sig" = "grey")) +
  labs(x = "H3K4me3 enrichment direction",
       y = "RNA-seq log2FC (NC vs FC)",
       title = "H3K4me3 enrichment correlates with gene expression",
       subtitle = "Peaks within 1500bp of TSS | OR=3.3, p=1.3e-11") +
  theme_bw()

integrated %>%
  filter(abs(distanceToTSS) < 1500, 
         H3K4_direction != "not_sig",
         RNA_direction != "not_DE") %>%
  count(H3K4_direction, RNA_direction) %>%
  mutate(concordant = case_when(
    H3K4_direction == "NC_enriched" & RNA_direction == "NC_up" ~ "Concordant",
    H3K4_direction == "FC_enriched" & RNA_direction == "FC_up" ~ "Concordant",
    TRUE ~ "Discordant"
  )) %>%
  ggplot(aes(x = H3K4_direction, y = n, fill = concordant)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Concordant" = "#2166ac", "Discordant" = "#d73027")) +
  labs(x = "H3K4me3 enrichment", 
       y = "Number of peaks",
       fill = "",
       title = "H3K4me3 concordance with RNA-seq",
       subtitle = "OR = 3.3, p = 1.3e-11 (Fisher's exact test)") +
  theme_bw()

colnames(integrated)
