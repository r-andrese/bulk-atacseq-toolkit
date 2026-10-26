## =============================================================================
## 03_differential_accessibility_and_annotation.R
##
## GC-aware normalisation (EDASeq), differential accessibility (edgeR) and
## peak annotation (ChIPseeker) of the consensus peak set from script 02.
##
## Differential testing: negative binomial GLM with a common dispersion and
## likelihood-ratio test (glmFit / glmLRT), test group vs reference group.
## =============================================================================

suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(EDASeq)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(txdbmaker)
  library(ChIPseeker)
  library(org.Mm.eg.db)
  library(ggupset)
})

## ---- editable parameters ----------------------------------------------------
base_dir      <- "/path/to/project"   # same as BASE in the shell scripts
consensus_dir <- file.path(base_dir, "results_atacseq/bwa/merged_library/macs2/broad_peak/consensus_by_group")
counts_file   <- file.path(consensus_dir, "counts_table/consensus_peaks.counts.tsv")
peaks_bed     <- file.path(consensus_dir, "consensus_peaks.bed")
fasta         <- file.path(base_dir, "references/Mus_musculus.GRCm38.dna_sm.toplevel.fa")  # needs .fai
gtf           <- file.path(base_dir, "references/Mus_musculus.GRCm38.95.gtf")   # same release as alignment
organism      <- "Mus musculus"
orgdb         <- "org.Mm.eg.db"
out_dir       <- file.path(base_dir, "results/consensus")

ref_group   <- "WT"     # sample names must start with "<group>_"
test_group  <- "HET"
bam_suffix  <- ".mLb.clN.sorted.bam"
chromosomes <- c(as.character(1:19), "X", "Y")

gc_bins     <- 20
fdr_cutoff  <- 0.05
tss_region  <- c(-3000, 3000)
## -----------------------------------------------------------------------------

qc_dir    <- file.path(out_dir, "QC")
annot_dir <- file.path(out_dir, "peak_annotation")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)


######################################################################
##  Differential accessibility
######################################################################

## ---- counts ----
cnt_table <- read.delim(counts_file, check.names = FALSE)
rownames(cnt_table) <- cnt_table$Geneid
count_cols <- 7:ncol(cnt_table)
colnames(cnt_table)[count_cols] <- sub(bam_suffix, "", basename(colnames(cnt_table)[count_cols]),
                                       fixed = TRUE)

cnt_table  <- cnt_table[cnt_table$Chr %in% chromosomes, ]
reads.peak <- as.matrix(cnt_table[, count_cols])

## ---- groups, assigned from sample names ----
samples <- colnames(reads.peak)
groups  <- ifelse(startsWith(samples, paste0(test_group, "_")), test_group,
           ifelse(startsWith(samples, paste0(ref_group,  "_")), ref_group, NA))
if (anyNA(groups))
  stop("Cannot assign a group to: ", paste(samples[is.na(groups)], collapse = ", "))
groups <- factor(groups, levels = c(ref_group, test_group))
print(table(groups))

## ---- GC content of peaks ----
gr <- GRanges(seqnames = cnt_table$Chr,
              ranges   = IRanges(cnt_table$Start, cnt_table$End))   # featureCounts: 1-based
gcContentPeaks <- letterFrequency(getSeq(FaFile(fasta), gr), "GC", as.prob = TRUE)[, 1]

## QC: log counts vs GC content per sample (lowess fit)
dfAll <- do.call(rbind, lapply(seq_len(ncol(reads.peak)), function(k) {
  fit <- lowess(x = gcContentPeaks, y = log1p(reads.peak[, k]), f = 1/10)
  data.frame(x = fit$x, y = fit$y, sample = samples[k])
}))
pdf(file.path(qc_dir, "GCcontent_peaks.pdf"))
print(ggplot(dfAll, aes(x, y, colour = sample)) + geom_line(linewidth = 1) +
        xlab("GC content") + ylab("log(count + 1)") + theme_classic())
dev.off()

## ---- GC-aware offsets (EDASeq) ----
eda <- newSeqExpressionSet(
  counts      = reads.peak,
  featureData = AnnotatedDataFrame(
    data.frame(gc = gcContentPeaks, row.names = rownames(reads.peak))))
eda <- withinLaneNormalization(eda, "gc", num.bins = gc_bins, which = "full", offset = TRUE)
eda <- betweenLaneNormalization(eda, which = "full", offset = TRUE)
dataOffset <- offst(eda)

## ---- edgeR with EDASeq offsets ----
design <- model.matrix(~groups)

d    <- DGEList(counts = reads.peak, group = groups)
keep <- filterByExpr(d)
d    <- d[keep, , keep.lib.sizes = FALSE]

d$offset   <- -dataOffset[keep, ]     # sign convention as in the EDASeq vignette
d          <- estimateGLMCommonDisp(d, design = design)
fit        <- glmFit(d, design = design)
lrt.EDASeq <- glmLRT(fit, coef = 2)   # test_group vs ref_group

DA_res <- as.data.frame(topTags(lrt.EDASeq, n = nrow(lrt.EDASeq$table)))
DA_res$Geneid <- rownames(DA_res)

## coordinates and per-group median raw counts
count_df <- data.frame(Geneid = rownames(d$counts))
count_df[[paste0("median_", ref_group,  "_counts")]] <- apply(d$counts[, groups == ref_group],  1, median)
count_df[[paste0("median_", test_group, "_counts")]] <- apply(d$counts[, groups == test_group], 1, median)

DA.res.coords <- DA_res %>%
  left_join(cnt_table[, c("Geneid", "Chr", "Start", "End")], by = "Geneid") %>%
  left_join(count_df, by = "Geneid")

write.table(DA.res.coords,
            file.path(out_dir, paste0("DA_", test_group, "_vs_", ref_group, ".tsv")),
            quote = FALSE, sep = "\t", row.names = FALSE)


######################################################################
##  Peak annotation
######################################################################

## TxDb built from the same Ensembl GTF used for alignment and quantification
txdb <- makeTxDbFromGFF(gtf, format = "gtf", organism = organism)

## helper: annotate a GRanges and write outputs
annotate_peaks <- function(gr, out) {
  if (length(gr) == 0) {
    message("No peaks for ", basename(out), "; annotation skipped")
    return(invisible(NULL))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ann <- annotatePeak(gr, tssRegion = tss_region, TxDb = txdb, annoDb = orgdb)

  df <- as.data.frame(ann)
  df$geneChr <- NULL
  colnames(df)[1:5] <- c("Chr", "peak_start", "peak_end", "peak_width", "peak_strand")
  write.table(df, file.path(out, "peaks_annotated.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(ann@annoStat, file.path(out, "annotation_summary.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  pdf(file.path(out, "AnnotVis.pdf")); print(upsetplot(ann, vennpie = TRUE)); dev.off()
  pdf(file.path(out, "TSSdist.pdf"), 6, 5)
  print(plotDistToTSS(ann, title = "ATAC-seq peaks relative to TSS")); dev.off()
  invisible(ann)
}

## all consensus peaks (BED is 0-based: start + 1 for GRanges)
peaks.bed <- read.table(peaks_bed, sep = "\t", header = FALSE)
peaks.gr  <- GRanges(seqnames = peaks.bed[, 1],
                     ranges   = IRanges(peaks.bed[, 2] + 1, peaks.bed[, 3]),
                     peakID   = peaks.bed[, 4])

pdf(file.path(annot_dir, "PeakCoverage.pdf"), 8, 20); print(covplot(peaks.gr)); dev.off()
annotate_peaks(peaks.gr, file.path(annot_dir, "all_consensus_peaks"))

## differentially accessible peaks
to_gr <- function(x) GRanges(seqnames = x$Chr, ranges = IRanges(x$Start, x$End),
                             peakID = x$Geneid)

annotate_peaks(to_gr(filter(DA.res.coords, FDR < fdr_cutoff, logFC > 0)),
               file.path(annot_dir, paste0("more_accessible_in_", test_group)))
annotate_peaks(to_gr(filter(DA.res.coords, FDR < fdr_cutoff, logFC < 0)),
               file.path(annot_dir, paste0("less_accessible_in_", test_group)))

## software versions
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
