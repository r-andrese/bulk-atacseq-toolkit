#!/bin/bash
# =============================================================================
# 01_run_nfcore_atacseq.sh
#
# Read QC, adapter trimming, BWA alignment, filtering and MACS2 broad peak
# calling with nf-core/atacseq v2.0.
#
# Notes
#   - The Nextflow head job needs adequate memory (32 GB was used).
# =============================================================================
set -euo pipefail

## ---- editable parameters ---------------------------------------------------
BASE=/path/to/project
REPO=/path/to/this/repository
REF=${BASE}/references

PIPELINE_VERSION=2.0
PROFILE=singularity
CONFIG=${REPO}/config/sge_cluster.config
SAMPLESHEET=${BASE}/download/samplesheet.csv
OUTDIR=${BASE}/results_atacseq
WORKDIR=${BASE}/work

FASTA=${REF}/Mus_musculus.GRCm38.dna_sm.toplevel.fa
GTF=${REF}/Mus_musculus.GRCm38.95.gtf
GENE_BED=${REF}/Mus_musculus.GRCm38.95.bed
BLACKLIST=${REF}/mm10-blacklist.GRCm38.bed
BWA_INDEX=${REF}/bwa

MITO_NAME=MT
MACS_GSIZE=2600000000
MIN_REPS_CONSENSUS=3   # nf-core's own consensus and DESeq2 QC only; not used downstream

EMAIL=""        # optional: address notified if the run fails
## ----------------------------------------------------------------------------

export TMPDIR=${BASE}/tmp
export NXF_TEMP=${TMPDIR}
export NXF_SINGULARITY_CACHEDIR=${BASE}/singularity_cache
mkdir -p "${TMPDIR}" "${NXF_SINGULARITY_CACHEDIR}"

EMAIL_OPT=()
if [ -n "${EMAIL}" ]; then EMAIL_OPT=(--email_on_fail "${EMAIL}"); fi

cd "${BASE}"

nextflow run nf-core/atacseq -r "${PIPELINE_VERSION}" \
  -profile "${PROFILE}" \
  -c "${CONFIG}" \
  --input "${SAMPLESHEET}" \
  --fasta "${FASTA}" \
  --gtf "${GTF}" \
  --bwa_index "${BWA_INDEX}" \
  --gene_bed "${GENE_BED}" \
  --blacklist "${BLACKLIST}" \
  --macs_gsize "${MACS_GSIZE}" \
  --mito_name "${MITO_NAME}" \
  --deseq2_vst \
  --min_reps_consensus "${MIN_REPS_CONSENSUS}" \
  --outdir "${OUTDIR}" \
  -work-dir "${WORKDIR}" \
  ${EMAIL_OPT[@]+"${EMAIL_OPT[@]}"} \
  -resume
