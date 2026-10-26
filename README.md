# bulk-atacseq-toolkit

Scripts for the analysis of bulk ATAC-seq data, from raw reads to differential
chromatin accessibility and peak annotation.

## Overview

This repository contains the workflow used to compare chromatin accessibility
between two experimental groups in bulk ATAC-seq data. Reads are processed and
peaks are called with the nf-core/atacseq pipeline (v2.0).

The downstream analysis proceeds as follows:

  - Consensus peak definition — Replicate-level MACS2 broad peaks are merged to
    generate a consensus peak set across samples.

  - Read counting — ATAC-seq fragments are quantified within each consensus peak.

  - Differential accessibility testing — Peak-level counts are normalised for GC
    content using EDASeq, then tested for differential accessibility between groups
    using edgeR.

  - Functional annotation — Differentially accessible peaks are annotated relative
    to gene models using ChIPseeker.

  - Motif enrichment — Transcription factor motif enrichment within differentially
    accessible peaks is assessed using HOMER.

The workflow is designed to identify genomic regions with altered chromatin
accessibility between conditions and to provide biological context through gene
annotation and regulatory motif analysis.

The scripts are provided as used for the analysis of mouse skin (heterozygous
vs wild-type; GEO accession GSE277902, mm10/GRCm38) reported in Veland et al.,
Nature (2026), doi:10.1038/s41586-026-11076-4.

- **Purpose:** differential chromatin accessibility between two groups of samples
- **Example dataset:** mouse skin, HET (n = 4) vs WT (n = 3), paired-end; GEO GSE277902
- **Genome and annotation:** GRCm38 (mm10), Ensembl release 95

## Repository structure

```
├── README.md
├── LICENSE
├── config/
│   ├── sge_cluster.config         Nextflow configuration (example, SGE cluster)
│   └── samplesheet_template.csv   nf-core/atacseq sample sheet with dummy values
└── scripts/
    ├── 00_prepare_inputs.sh
    ├── 01_run_nfcore_atacseq.sh
    ├── 02_consensus_peaks_featurecounts.sh
    ├── 03_differential_accessibility_and_annotation.R
    ├── 04_motif_enrichment_homer.sh
    └── 05_plot_motif_enrichment.R
```

## Workflow

1. **`00_prepare_inputs.sh`** converts raw reads from SRA to
   gzipped FASTQ, generates the nf-core sample sheet from the SRA and GEO
   metadata, and prepares the reference files.
   - Sample titles must follow `<group>_REP<n>`; several SRA runs of the same
     sample are entered as technical replicates (same `sample` and `replicate`)
2. **`01_run_nfcore_atacseq.sh`** runs nf-core/atacseq v2.0: read QC and
   trimming, alignment with BWA-MEM, filtering, merging of technical
   replicates, and MACS2 broad peak calling.
3. **`02_consensus_peaks_featurecounts.sh`** defines consensus peaks per group
   and counts read pairs per consensus peak with featureCounts.
   - A merged region is retained if it is overlapped, with at least 50 %
     reciprocal overlap, by peaks from ≥ 2 of 3 WT or ≥ 3 of 4 HET replicates
   - The two group-level sets are combined.
4. **`03_differential_accessibility_and_annotation.R`** normalises counts,
   tests for differential accessibility and annotates peaks.
   - Within-sample GC-content normalisation followed by between-sample
     normalisation (EDASeq, full quantile), applied as offsets in edgeR
   - Negative binomial GLM with common dispersion, likelihood-ratio test
   - Annotation of all consensus peaks and of significant peaks (FDR < 0.05),
     separately for gain and loss of accessibility
5. **`04_motif_enrichment_homer.sh`** tests known and de novo motif
   enrichment (HOMER) in peaks gaining and losing accessibility.
6. **`05_plot_motif_enrichment.R`** compares known-motif enrichment between
   gained and lost peaks.

## Inputs

- Paired-end ATAC-seq FASTQ files (or SRA accessions)
- `SraRunTable.csv` from the SRA Run Selector and the GEO series matrix
  (only needed by `00_prepare_inputs.sh`)
- Genome FASTA (indexed with `samtools faidx`) and GTF from the same Ensembl release
- ENCODE blacklist for the genome assembly
- Sample sheet in the nf-core/atacseq v2.0 format (`sample,fastq_1,fastq_2,replicate`)

## Outputs

Paths are relative to the project directory (`BASE`).

- `results_atacseq/`: nf-core/atacseq output (alignments, MACS2 peaks, MultiQC)
- `results_atacseq/bwa/merged_library/macs2/broad_peak/consensus_by_group/`:
  consensus peaks (`consensus_peaks.bed`, `consensus_peaks.saf`) and the count
  matrix (`counts_table/consensus_peaks.counts.tsv`)

In `results/consensus/`:

- `DA_<test>_vs_<reference>.tsv`: edgeR results with peak coordinates and
  per-group median raw counts
- `QC/GCcontent_peaks.pdf`: per-sample relationship between counts and GC content
- `peak_annotation/`: annotated peak tables, annotation summaries and plots
  for all consensus peaks and for peaks gaining or losing accessibility
- `sessionInfo.txt`: R and package versions
- `motif_enrichment/`: HOMER results for gained and lost peaks, and the
  gained-vs-lost enrichment plot and table


## Requirements

- Nextflow and Singularity (or another container engine supported by nf-core)
- SRA Toolkit, wget, GNU awk, samtools, BWA
- bedtools, BEDOPS, Subread (featureCounts, version 2.0.2 or later)
- HOMER
- R with Bioconductor: EDASeq, edgeR, ChIPseeker, txdbmaker, org.Mm.eg.db,
  Rsamtools, Biostrings, GenomicRanges; CRAN: dplyr, ggplot2, ggupset, ggrepel

## Usage

1. Edit the parameter block at the top of each script.
2. Run the scripts in order:

```bash
bash scripts/00_prepare_inputs.sh
bash scripts/01_run_nfcore_atacseq.sh          # submit as a long-running job
bash scripts/02_consensus_peaks_featurecounts.sh
Rscript scripts/03_differential_accessibility_and_annotation.R
bash scripts/04_motif_enrichment_homer.sh
Rscript scripts/05_plot_motif_enrichment.R
```

## Data availability

Raw data for the example analysis are available from GEO under accession
[GSE277902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277902).

## Citations

- Veland N *et al.* γδ T cell receptor dependencies define a unique immunosurveillance
  modality. *Nature* (2026). doi:10.1038/s41586-026-11076-4
- Ewels PA *et al.* The nf-core framework for community-curated bioinformatics
  pipelines. *Nat Biotechnol* 38, 276–278 (2020).
- nf-core/atacseq: doi:[10.5281/zenodo.2634132](https://doi.org/10.5281/zenodo.2634132)
- Zhang Y *et al.* Model-based analysis of ChIP-Seq (MACS). *Genome Biol* 9, R137 (2008).
- Liao Y, Smyth GK, Shi W. featureCounts. *Bioinformatics* 30, 923–930 (2014).
- Risso D *et al.* GC-content normalization for RNA-Seq data. *BMC Bioinformatics* 12, 480 (2011).
- Robinson MD, McCarthy DJ, Smyth GK. edgeR. *Bioinformatics* 26, 139–140 (2010).
- Yu G, Wang LG, He QY. ChIPseeker. *Bioinformatics* 31, 2382–2383 (2015).
- Heinz S *et al.* Simple combinations of lineage-determining transcription factors
  prime cis-regulatory elements required for macrophage and B cell identities.
  *Mol Cell* 38, 576–589 (2010).
- Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist. *Sci Rep* 9, 9354 (2019).

## Licence

Released under the MIT Licence; see [LICENSE](LICENSE).

## Author

Rosa Andres Ejarque, King's College London
