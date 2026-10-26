#!/bin/bash
# =============================================================================
# 02_consensus_peaks_featurecounts.sh
#
# Builds a consensus peak set from the nf-core/atacseq MACS2 broad peaks output
# and counts read pairs per peak with featureCounts.
#
# Consensus definition, per group:
#   1. merge the peaks of all replicates in the group (bedtools merge);
#   2. keep merged regions overlapped, with >= MIN_OVERLAP reciprocal overlap,
#      by at least REF_MIN_REPS / TEST_MIN_REPS peaks across the group's
#      replicate files (equivalent to replicates, as MACS2 peaks within one
#      sample do not overlap);
#   3. take the union of the two group-level sets (bedops --merge).
#
# Outputs (in ${OUT}):
#   consensus_peaks.saf          featureCounts annotation (1-based, SAF)
#   consensus_peaks.bed          peaks with IDs (0-based, standard BED)
#   counts_table/consensus_peaks.counts.tsv
#
# Paths must not contain spaces (sample lists are split on whitespace).
#
# Requirements: bedtools, bedops, Subread featureCounts >= 2.0.2
#   (--countReadPairs counts fragments; with Subread < 2.0.2, -p alone did so)
# =============================================================================
set -euo pipefail
export LC_ALL=C

## ---- editable parameters ---------------------------------------------------
NFCORE_OUTDIR=/path/to/project/results_atacseq
ROOT=${NFCORE_OUTDIR}/bwa/merged_library
PEAKDIR=${ROOT}/macs2/broad_peak
OUT=${PEAKDIR}/consensus_by_group   # kept apart from nf-core's own consensus/ folder

PEAK_SUFFIX=.mLb.clN_peaks.broadPeak
BAM_SUFFIX=.mLb.clN.sorted.bam

REF_GROUP=WT
REF_SAMPLES="WT_Skin_REP1 WT_Skin_REP2 WT_Skin_REP3"
REF_MIN_REPS=2

TEST_GROUP=HET
TEST_SAMPLES="HET_Skin_REP1 HET_Skin_REP2 HET_Skin_REP3 HET_Skin_REP4"
TEST_MIN_REPS=3

MIN_OVERLAP=0.50              # reciprocal overlap fraction (bedtools -f -r)
PEAK_ID_PREFIX=merged_overlap_macs2_
FC_FRAC_OVERLAP=0.2           # featureCounts --fracOverlap
FC_STRAND=0                   # featureCounts -s (0 = unstranded)
THREADS=16
## ----------------------------------------------------------------------------

mkdir -p "${OUT}/counts_table"

peak_files() { for s in "$@"; do echo "${PEAKDIR}/${s}${PEAK_SUFFIX}"; done; }

for f in $(peak_files ${REF_SAMPLES} ${TEST_SAMPLES}); do
  [ -s "${f}" ] || { echo "Missing peak file: ${f}" >&2; exit 1; }
done

# $1 label, $2 minimum number of replicates, $3... sample names
group_consensus() {
  local label=$1 min=$2; shift 2
  local files; files=$(peak_files "$@")
  cat ${files} | sort -k1,1 -k2,2n | bedtools merge -i - > "${OUT}/merged_${label}.bed"
  bedtools intersect -a "${OUT}/merged_${label}.bed" -b ${files} \
      -f "${MIN_OVERLAP}" -r -c -nonamecheck \
    | awk -v min="${min}" '$NF >= min' > "${OUT}/consensus_${label}.bed"
}

group_consensus "${REF_GROUP}"  "${REF_MIN_REPS}"  ${REF_SAMPLES}
group_consensus "${TEST_GROUP}" "${TEST_MIN_REPS}" ${TEST_SAMPLES}

bedops --merge "${OUT}/consensus_${REF_GROUP}.bed" "${OUT}/consensus_${TEST_GROUP}.bed" \
  > "${OUT}/consensus_union.bed"

# SAF is 1-based (start + 1); the BED keeps the 0-based start. IDs match.
awk -F'\t' -v OFS='\t' -v p="${PEAK_ID_PREFIX}" \
  '{ print p NR, $1, $2 + 1, $3, "." }' "${OUT}/consensus_union.bed" > "${OUT}/consensus_peaks.saf"
awk -F'\t' -v OFS='\t' -v p="${PEAK_ID_PREFIX}" \
  '{ print $1, $2, $3, p NR, 0, "." }' "${OUT}/consensus_union.bed" > "${OUT}/consensus_peaks.bed"

BAMS=()
for s in ${TEST_SAMPLES} ${REF_SAMPLES}; do BAMS+=("${ROOT}/${s}${BAM_SUFFIX}"); done

featureCounts -F SAF -O --fracOverlap "${FC_FRAC_OVERLAP}" \
  -p --countReadPairs -s "${FC_STRAND}" -T "${THREADS}" \
  -a "${OUT}/consensus_peaks.saf" \
  -o "${OUT}/counts_table/consensus_peaks.counts" "${BAMS[@]}"

# drop featureCounts' leading comment line so R reads the table directly
tail -n +2 "${OUT}/counts_table/consensus_peaks.counts" \
  > "${OUT}/counts_table/consensus_peaks.counts.tsv"
