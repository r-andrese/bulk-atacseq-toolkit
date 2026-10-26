## =============================================================================
## 05_plot_motif_enrichment.R
##
## Compares known-motif enrichment in peaks gaining vs losing accessibility
## (HOMER knownResults.txt from script 04): one point per motif, with
## -log(FDR) in gained peaks on the x axis and in lost peaks on the y axis.
## FDR is Benjamini-Hochberg on HOMER's known-motif p-values, per peak set.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

## ---- editable parameters ----------------------------------------------------
base_dir   <- "/path/to/project"   # same as BASE in the shell scripts
motif_dir  <- file.path(base_dir, "results/consensus/motif_enrichment")
test_group <- "HET"
out_pdf    <- file.path(motif_dir, paste0("motif_enrichment_gained_vs_lost_", test_group, ".pdf"))

## Motifs to label, matched to the short name.
## NULL labels the n_top motifs with the strongest enrichment in either set.
highlight_motifs <- NULL
n_top            <- 15
## -----------------------------------------------------------------------------

read_known <- function(set) {
  f <- file.path(motif_dir, paste0("homer_", set, "_in_", test_group), "knownResults.txt")
  if (!file.exists(f))
    stop("Missing ", f, ": script 04 found no ", set, " peaks or did not run; ",
         "the comparison needs both peak sets.")
  df <- read.delim(f, check.names = TRUE) %>% select(Motif.Name, P.value)
  df$FDR <- p.adjust(df$P.value, method = "BH")
  # p-values below double precision are read as 0; floor them to keep -log finite
  df$logFDR <- -log(pmax(df$FDR, .Machine$double.xmin))
  df
}

ggdf <- merge(read_known("gained"), read_known("lost"),
              by = "Motif.Name", suffixes = c(".gained", ".lost"))

ggdf$motif <- sub("/.*", "", gsub("\\(.*?\\)", "", ggdf$Motif.Name))

if (is.null(highlight_motifs)) {
  top <- order(pmax(ggdf$logFDR.gained, ggdf$logFDR.lost), decreasing = TRUE)
  highlight_motifs <- ggdf$motif[head(top, n_top)]
}
ggdf$highlight <- ggdf$motif %in% highlight_motifs

pdf(out_pdf, 6, 6)
print(
  ggplot(ggdf, aes(logFDR.gained, logFDR.lost)) +
    geom_point(aes(colour = highlight)) +
    geom_text_repel(data = subset(ggdf, highlight), aes(label = motif),
                    colour = "black", box.padding = 0.6, max.overlaps = Inf) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    theme_bw() +
    labs(x = paste0("Enrichment in peaks gained in ", test_group, " (-ln FDR)"),
         y = paste0("Enrichment in peaks lost in ", test_group, " (-ln FDR)"))
)
dev.off()

write.table(ggdf, file.path(motif_dir, paste0("motif_enrichment_gained_vs_lost_", test_group, ".tsv")),
            quote = FALSE, sep = "\t", row.names = FALSE)
