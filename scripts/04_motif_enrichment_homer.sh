#!/bin/bash
# =============================================================================
# 04_motif_enrichment_homer.sh
#
# Known and de novo motif enrichment (HOMER findMotifsGenome.pl) in peaks that
# gain or lose accessibility, using the differential accessibility table from
# script 03.
#
# The genome is given as the FASTA used for alignment rather than a HOMER
# genome package.
# The FASTA is soft-masked, so -mask excludes repeats.
#
# Background (USE_CONSENSUS_BACKGROUND):
#   false (default)  HOMER's GC-matched random genomic background
#   true             all consensus peaks (script 02), i.e. enrichment among
#                    changed regions relative to accessible regions
#
# Requirements: HOMER (findMotifsGenome.pl), GNU awk (gawk)
# =============================================================================
set -euo pipefail

## ---- editable parameters ---------------------------------------------------
BASE=/path/to/project
DA_TABLE=${BASE}/results/consensus/DA_HET_vs_WT.tsv
CONSENSUS_BED=${BASE}/results_atacseq/bwa/merged_library/macs2/broad_peak/consensus_by_group/consensus_peaks.bed
FASTA=${BASE}/references/Mus_musculus.GRCm38.dna_sm.toplevel.fa
OUT=${BASE}/results/consensus/motif_enrichment
PREPARSED_DIR=${BASE}/homer_preparsed      # HOMER background cache

TEST_GROUP=HET
FDR_CUTOFF=0.05
SIZE=200                                   # bp around the peak centre
THREADS=16
USE_CONSENSUS_BACKGROUND=false
## ----------------------------------------------------------------------------

mkdir -p "${OUT}" "${PREPARSED_DIR}"

# Peaks gaining (logFC > 0) or losing (logFC < 0) accessibility, as BED.
# Columns are located by name; DA coordinates are 1-based, BED start is 0-based.
make_bed() {   # $1 = "gained" | "lost"
  gawk -F'\t' -v OFS='\t' -v dir="$1" -v fdr="${FDR_CUTOFF}" '
  NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
  $col["FDR"] < fdr && ((dir == "gained" && $col["logFC"] > 0) || (dir == "lost" && $col["logFC"] < 0)) {
    print $col["Chr"], $col["Start"] - 1, $col["End"], $col["Geneid"], 0, "+"
  }' "${DA_TABLE}"
}

BG_OPT=()
if [ "${USE_CONSENSUS_BACKGROUND}" = true ]; then
  gawk -v OFS='\t' '{ print $1, $2, $3, $4, 0, "+" }' "${CONSENSUS_BED}" > "${OUT}/background_consensus.bed"
  BG_OPT=(-bg "${OUT}/background_consensus.bed")
fi

for dir in gained lost; do
  bed="${OUT}/${dir}_in_${TEST_GROUP}.bed"
  make_bed "${dir}" > "${bed}"
  n=$(wc -l < "${bed}")
  echo "${dir}: ${n} peaks"
  if [ "${n}" -eq 0 ]; then echo "  no peaks, skipping"; continue; fi

  findMotifsGenome.pl "${bed}" "${FASTA}" "${OUT}/homer_${dir}_in_${TEST_GROUP}" \
    -size "${SIZE}" -mask -p "${THREADS}" -preparsedDir "${PREPARSED_DIR}" \
    ${BG_OPT[@]+"${BG_OPT[@]}"}
done
