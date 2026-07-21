setwd("~/Dropbox/Lab/Projects/Beetles/Genome_Reassembly/GRN_ananse/peaks/peaks_mapq13/")

library(DESeq2)
library(tidyverse)
library(janitor)

## We can annotate peaks and look at distance to genes using homer
### lets do some comparisons of Genrich combined vs Genrich condition-separate
## remove chrMT and scaffolds
## grep -v "^chrMT" larvae.mapq13.bed | grep -v "^scaffold_" > larvae.mapq13.filtered.bed
## Use HOMER to annotate peaks
## annotatePeaks.pl larvae.mapq13.bed nicrophorus > larvae_annotated.txt

larvae_annot <- read_tsv("larvae_annotated.txt")
larvae_annot <- clean_names(larvae_annot)

head(larvae_annot)

larvae_annot <- larvae_annot %>%
  rename(PeakID = 1) %>%
  mutate(Feature = case_when(
    grepl("promoter", annotation, ignore.case = TRUE) ~ "Promoter",
    grepl("exon", annotation, ignore.case = TRUE)     ~ "Exon",
    grepl("intron", annotation, ignore.case = TRUE)   ~ "Intron",
    grepl("TTS", annotation, ignore.case = TRUE)      ~ "TTS",
    grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Intergenic",
    TRUE ~ "Other"
  ))

# Peak counts per feature
larvae_peak_counts <- larvae_annot %>%
  count(Feature, name = "n_peaks")

# Genome feature sizes from HOMER's console output
genome_sizes <- tibble(
  Feature      = c("Promoter", "Exon", "Intron", "Intergenic", "TTS"),
  total_bp     = c(15364045, 17146615, 58463858, 99375058, 11237606)
)

genome_size_total <- sum(genome_sizes$total_bp)
larvae_total_peaks <- nrow(larvae_annot)

larvae_summary_table <- larvae_peak_counts %>%
  left_join(genome_sizes, by = "Feature") %>%
  mutate(
    expected   = (total_bp / genome_size_total) * larvae_total_peaks,
    pct_peaks  = (n_peaks / larvae_total_peaks) * 100,
    log2_ratio = log2(n_peaks / expected),
    logP       = pbinom(n_peaks, larvae_total_peaks, total_bp / genome_size_total, 
                        lower.tail = FALSE, log.p = TRUE) / log(10)
  )

# LogP enrichment (+values depleted)
larvae_summary_table


#### DEseq2

## Here the counts have been donne in one go with all samples in genrich
## mapq13 filtered, q=0.05

counts <- read.table("larvae.mapq13.counts") %>% 
  select(-c(V5,V6))

head(counts)

colnames(counts) <- c("seqnames",
                      "start", 
                      "end",
                      "peak",
                      "A3_FC", "A4_FC", "A5_FC", "A6_FC", 
                      "A3_NC", "A4_NC", "A5_NC", "A6_NC")

head(counts)


## remove chrMT and unplaced scaffolds
keep_chr <- !grepl("^(chrMT|scaffold_)", counts$seqnames)
counts <- counts[keep_chr, ]

## start DESeq2
countData <- counts %>% 
  select(4:12) %>% 
  column_to_rownames(var="peak")


countData <- as.matrix(countData)

metaData = data.frame(ID = c("A3_FC", "A4_FC", "A5_FC", "A6_FC", 
                             "A3_NC", "A4_NC", "A5_NC", "A6_NC"))
metaData$Cond = as.factor(substr(metaData$ID,4,5))
metaData$Gen = as.factor(rep(c("g1", "g2", "g3", "g4")))

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~Gen + Cond)

#remove features with low counts
dds = dds[rowMeans(counts(dds)) > 10, ]
nrow(dds) # 27149

# rlog transform counts (by average length and correcting for library size)
rld = rlog(dds, blind=FALSE)

plotPCA(rld, intgroup = "Cond", ntop = nrow(rld))

plotPCA(rld, intgroup = "Cond", ntop = nrow(rld)) +
  geom_text(aes(label = rld$ID), nudge_y = 1, size = 3)

plotPCA(rld, intgroup = "Cond", ntop = 2000) +
  theme_bw()

plotPCA(rld, intgroup = "Cond", ntop = 2000) +
  geom_text(aes(label = rld$ID), nudge_y = 1, size = 3)


vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Cond", ntop = 10000) +
  geom_text(aes(label = rld$ID), nudge_y = 1, size = 3)


colSums(counts(dds))

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Cond_NC_vs_FC")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="Cond_NC_vs_FC", type="apeglm")

length(which(res$padj < 0.05))
#315

length(which(res$padj < 0.10))
#465

peaks.dds <- data.frame(res) %>% 
  rownames_to_column(var = "peak")

head(peaks.dds)

all.peaks <- counts %>% 
  dplyr::select(seqnames, start, end, peak) %>% 
  filter(peak %in% peaks.dds$peak) %>% 
  left_join(peaks.dds)

write_tsv(all.peaks, "DA_res_allpeaks.tsv")

head(all.peaks)

head(counts)

sig.dds <- data.frame(res) %>% 
  rownames_to_column(var = "peak") %>% 
  filter(padj < 0.05)

sig.peaks <- counts %>% 
  dplyr::select(seqnames, start, end, peak) %>% 
  filter(peak %in% sig.dds$peak) %>% 
  left_join(peaks.dds)

## make an all sig peaks file to annotate in homer
sig.peaks %>% 
  select(seqnames,start,end,peak) %>% 
  write_tsv(., "larvae_differentially_acccesible.bed", col_names = F)

  

## make a file for accessible in NC 
sig.peaks %>% 
  filter(log2FoldChange > 0) %>% 
  dplyr::select(seqnames,start,end,peak) %>% 
  write_tsv(., "~/homer_res/NC_acccesible.bed", col_names = F)

sig.peaks %>% 
  filter(log2FoldChange < 0) %>% 
  dplyr::select(seqnames,start,end,peak) %>% 
  write_tsv(., "~/homer_res/FC_acccesible.bed", col_names = F)

all.peaks %>% 
  drop_na(padj) %>% 
  dplyr::select(seqnames,start,end,peak) %>% 
  write_tsv(., "~/homer_res/bg_peaks.bed", col_names = F)
  

####### PLOTS #####

## get differentially accessible annotations
## annotatePeaks.pl larvae_differentially_acccesible.bed nicrophorus > larvae_differentially_accessible_annotated.txt

# Load DA peaks annotation
da_annot <- read_tsv("larvae_differentially_accessible_annotated.txt")
da_annot <- clean_names(da_annot)
da_annot <- da_annot %>%
  rename(PeakID = 1) %>%
  mutate(Feature = case_when(
    grepl("promoter", annotation, ignore.case = TRUE) ~ "Promoter",
    grepl("exon", annotation, ignore.case = TRUE)     ~ "Exon",
    grepl("intron", annotation, ignore.case = TRUE)   ~ "Intron",
    grepl("TTS", annotation, ignore.case = TRUE)      ~ "TTS",
    grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Intergenic",
    TRUE ~ "Other"
  ))

# Combine with label
combined_annot <- bind_rows(
  larvae_annot %>% mutate(set = "All OCRs"),
    da_annot %>% mutate(set = "DA OCRs")
)

## plot cumulative density
ggplot(larvae_annot, aes(x = distance_to_tss)) +
  stat_ecdf() +
  xlim(-50000, 50000) +
  labs(x = "Distance to gene (bp)", y = "Cumulative density",
       title = "Cumulative density of OCR distance to TSS") +
  theme_classic()

larvae_annot %>%
  filter(abs(distance_to_tss) <= 50000) %>%
  ggplot(aes(x = distance_to_tss)) +
  stat_ecdf() +
  geom_vline(xintercept = c(-30000, -15000, -6000, -1000, 1000, 6000, 15000, 30000), colour = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = c(0,1), linetype = "dashed", colour = "grey60", alpha = 0.7) +
  geom_vline(xintercept = c(0), linetype = "dashed", colour = "grey60", alpha = 0.7) +
  labs(x = "Distance to gene (bp)", y = "Cumulative density",
       title = "Cumulative density of OCR distance to TSS") +
  xlim(-45000,45000) +
  theme_classic()

########

# Feature plot
feat_plot <- combined_annot %>%
  filter(Feature != "Other") %>%
  count(set, Feature) %>%
  group_by(set) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = set, y = pct, fill = Feature)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "", y = "Fraction of OCRs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Distance bins plot
bin_plot <- combined_annot %>%
  mutate(bin = cut(abs(distance_to_tss),
                   breaks = c(0, 1000, 6000, 15000, 30000, Inf),
                   labels = c("0-1Kb", "1-6Kb", "6Kb-15Kb", "15-30Kb", "beyond 30Kb"),
                   include.lowest = TRUE)) %>%
  count(set, bin) %>%
  group_by(set) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = set, y = pct, fill = bin)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Greens") +
  labs(x = "", y = "Fraction of OCRs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

bin_plot + feat_plot

# library(RColorBrewer)
# # See all 9 shades in the palette
# brewer.pal(9, "Blues")
# # Take shades 2-6 (skipping the lightest)
# my_purples <- brewer.pal(9, "Blues")[2:6]
# # Then use in your plot

### peaks per gene

raw_peaks_per_gene <- larvae_annot %>%
  count(nearest_promoter_id, name = "n_peaks")

head(raw_peaks_per_gene)
peaks_per_gene <- larvae_annot %>%
  count(nearest_promoter_id, name = "n_peaks")

head(peaks_per_gene)

# get all gene IDs from GTF
# grep -v "^#" ~/new_hic_genome/*.gtf | \
# awk '$3=="mRNA"' | \
# grep -o 'transcript_id "[^"]*"' | \
# sed 's/transcript_id "//;s/"//' | \
# sort -u > all_genes.txt

all_genes <- read.table("all_genes.txt", col.names = "nearest_promoter_id")

peaks_per_gene <- all_genes %>%
  left_join(raw_peaks_per_gene, by = "nearest_promoter_id") %>%
  replace_na(list(n_peaks = 0)) %>%
  mutate(bin = case_when(
    n_peaks == 0  ~ "0",
    n_peaks == 1  ~ "1",
    n_peaks == 2  ~ "2",
    n_peaks <= 5  ~ "3-5",
    n_peaks <= 10 ~ "6-10",
    TRUE          ~ ">10"
  ))

table(peaks_per_gene$bin)

table(peaks_per_gene$bin)

bin_levels <- c("0", "1", "2", "3-5", "6-10", ">10")

peaks_per_gene %>%
  count(bin) %>%
  mutate(
    pct = n / sum(n),
    bin = factor(bin, levels = bin_levels)
  ) %>%
  ggplot(aes(x = "Larvae", y = pct, fill = bin)) +
  geom_col() +
  scale_fill_manual(values = c("white", "#e5f5e0", "#a1d99b", "#41ae76", "#238b45", "#00441b"),
                    breaks = rev(bin_levels)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Fraction of genes", 
       title = "Fraction of genes by number of associated OCRs",
       fill = "Number of OCRs") +
  theme_classic()


# DA peaks per gene
raw_da_peaks_per_gene <- da_annot %>%
  count(nearest_promoter_id, name = "n_peaks")

# Build both with labels
all_per_gene <- all_genes %>%
  left_join(raw_peaks_per_gene, by = "nearest_promoter_id") %>%
  replace_na(list(n_peaks = 0)) %>%
  mutate(set = "All OCRs")

da_per_gene <- all_genes %>%
  left_join(raw_da_peaks_per_gene, by = "nearest_promoter_id") %>%
  replace_na(list(n_peaks = 0)) %>%
  mutate(set = "Differentially Accessible OCRs")

# Combine and bin
combined_per_gene <- bind_rows(all_per_gene, da_per_gene) %>%
  mutate(bin = case_when(
    n_peaks == 0  ~ "0",
    n_peaks == 1  ~ "1",
    n_peaks == 2  ~ "2",
    n_peaks <= 5  ~ "3-5",
    n_peaks <= 10 ~ "6-10",
    TRUE          ~ ">10"
  ))

table(combined_per_gene$set, combined_per_gene$bin)

# Plot
bin_levels <- c("0", "1", "2", "3-5", "6-10", ">10")

combined_per_gene %>%
  filter(n_peaks > 0) %>%
  count(set, bin) %>%
  group_by(set) %>%
  mutate(pct = n / sum(n),
         bin = factor(bin, levels = bin_levels)) %>%
  ggplot(aes(x = set, y = pct, fill = bin)) +
  geom_col() +
  scale_fill_manual(values = c("white", "#e5f5e0", "#a1d99b", "#41ae76", "#238b45", "#00441b"),
                    breaks = rev(bin_levels)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Fraction of genes",
       title = "Fraction of genes by number of associated OCRs",
       fill = "Number of OCRs") +
  theme_classic()




summary(abs(da_annot$distance_to_tss))
hist(abs(da_annot$distance_to_tss), breaks = 50)


hist(abs(larvae_annot$distance_to_tss), breaks = 50)

## rna-seq results
rna_de <- read_tsv("../../../NewRNAseq/results/Care_Deseq2.tsv")
head(rna_de)

rna_de <- read_tsv("../../../NewRNAseq/results/Care_Deseq2.tsv") %>% 
  mutate(isDE = ifelse(padj < 0.05, 1, 0)) %>% 
  select(gene_id, gene_name, isDE)

da_annot_full <- larvae_annot %>%
  filter(abs(distance_to_tss) <= 20000) %>%
  inner_join(all.peaks, by = c("PeakID" = "peak")) %>%
  mutate(gene_id = sub("\\.t\\d+$", "", nearest_promoter_id)) %>% 
  mutate(isDA = ifelse(padj < 0.05, 1, 0)) %>% 
  select(PeakID,chr,Feature,distance_to_tss,gene_id,isDA) %>% 
  left_join(rna_de)

da_annot_full

fisher.test(da_annot_full$isDE, da_annot_full$isDA)

table(da_annot_full$isDA, da_annot_full$isDE)

table(da_annot_full$isDA, da_annot_full$isDE)

da_annot_full %>%
  group_by(Feature) %>%
  summarise(
    total = n(),
    n_DA = sum(isDA, na.rm = TRUE),
    n_DE = sum(isDE, na.rm = TRUE),
    n_DA_and_DE = sum(isDA == 1 & isDE == 1, na.rm = TRUE)
  )

da_annot_full %>%
  filter(!is.na(isDE), !is.na(isDA)) %>%
  group_by(Feature) %>%
  summarise(
    fisher_p = fisher.test(isDA, isDE)$p.value,
    odds_ratio = fisher.test(isDA, isDE)$estimate
  ) %>%
  arrange(fisher_p)
## make a background file of all peaks

### annotate peaks ###
### tihis is just playing around trying to identify 'shadow enhancers' ###
## code is a bit of a mess ###

# 1. Peaks with DE results
peaks.dds <- data.frame(res) %>% 
  rownames_to_column(var = "peak")

all.peaks <- counts %>% 
  dplyr::select(seqnames, start, end, peak) %>% 
  filter(peak %in% peaks.dds$peak) %>% 
  left_join(peaks.dds, by = "peak")

# 2. Make GRanges from peaks
peak_gr <- GRanges(
  seqnames = all.peaks$seqnames,
  ranges = IRanges(start = all.peaks$start, end = all.peaks$end),
  peak = all.peaks$peak
)

# 3. Get gene TSSs
txdb <- makeTxDbFromGFF("~/new_hic_genome/Nvi2.liftoff.mitohifi.final.ok.genename.for.cellranger.gtf")
g <- genes(txdb, single.strand.genes.only = FALSE)
if (is(g, "GRangesList")) {
  g <- unlist(range(g))
}
tss <- promoters(g, upstream = 0, downstream = 1)

# 4. Find nearest TSS for each peak
hits <- distanceToNearest(peak_gr, tss)

# 5. Add gene annotation to peaks
all.peaks$geneId <- NA
all.peaks$distanceToTSS <- NA
all.peaks$geneId[queryHits(hits)] <- names(tss)[subjectHits(hits)]
all.peaks$distanceToTSS[queryHits(hits)] <- mcols(hits)$distance

# 6. Filter for distal peaks
distal <- all.peaks %>%
  filter(distanceToTSS > 2000)

# 7. Shadow enhancer identification
enhancer_counts <- table(distal$gene_name)
shadow_genes <- names(enhancer_counts[enhancer_counts >= 2])

shadow_peaks <- distal %>%
  filter(gene_name %in% shadow_genes)

shadow_summary <- shadow_peaks %>%
  group_by(gene_name) %>%
  summarise(
    n_enhancers = n(),
    n_sig = sum(padj < 0.10, na.rm = TRUE),
    n_stable = sum(padj > 0.10 | is.na(padj)),
    mean_lfc = mean(log2FoldChange, na.rm = TRUE)
  )

buffered <- shadow_summary %>%
  filter(n_stable >= 1 & n_sig >= 1)

nrow(buffered)
buffered

DiffExp %>% filter(gene_name == "MYC")
fc_enh[target == "MYC"]
nc_enh[target== "MYC"]

# Read enhancer-only networks
fc_enh <- fread("../../results2/fc_enhancer_influence_diffnetwork.tsv")
nc_enh <- fread("../../results2/NC_enhancer_influence_diffnetwork.tsv")

head(nc_enh)

# Filter forbuffered shadow enhancer genes
buffered_genes <- buffered$gene_name
fc_shadow <- fc_enh[target %in% buffered_genes]
nc_shadow <- nc_enh[target %in% buffered_genes]

## find number of TFs regulating buffered genes in the diff network
fc_shadow_tfs <- fc_shadow %>% distinct(source) %>% pull
nc_shadow_tfs <- nc_shadow %>% distinct(source) %>% pull

head(fc_enh)
head(nc_enh)
colnames(fc_enh)
nrow(fc_shadow)
nrow(nc_shadow)

length(fc_shadow_tfs)
length(nc_shadow_tfs)
# Overlap?
length(intersect(fc_shadow_tfs, nc_shadow_tfs))


# Total genes in diffnetwork
all_fc_targets <- fc_enh %>% distinct(target) %>% pull
all_nc_targets <- nc_enh %>% distinct(target) %>% pull

# Proportion of shadow genes in diffnetwork vs background
fisher.test(matrix(c(
  length(unique(fc_shadow$target)),
  length(buffered_genes) - length(unique(fc_shadow$target)),
  length(all_fc_targets) - length(unique(fc_shadow$target)),
  nrow(shadow_summary) - length(all_fc_targets)
), nrow = 2))


# TFs unique to each condition
tf_fc_unique <- setdiff(fc_shadow_tfs, nc_shadow_tfs)
tf_nc_unique <- setdiff(nc_shadow_tfs, fc_shadow_tfs)


nc_enh_influence <- read_tsv("../../results2/nc_enhancer_influence.tsv", show_col_types = FALSE)
fc_enh_influence <- read_tsv("../../results2/fc_enhancer_influence.tsv", show_col_types = FALSE)


nc_enh_influence %>% 
  filter(factor %in% tf_nc_unique ) %>%
  dplyr::select(factor, influence_score, direct_targets, factor_fc) %>%
  print(n=27)

fc_enh_influence %>% 
  filter(factor %in% tf_fc_unique ) %>%
  dplyr::select(factor, influence_score, direct_targets, factor_fc) %>%
  print(n=27)

nc_enh_influence %>% filter(factor == "HSF")
fc_enh_influence %>% filter(factor == "HSF")

shadow_peaks %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up_in_NC", "down_in_NC")) %>%
  pull(direction) %>%
  table()

