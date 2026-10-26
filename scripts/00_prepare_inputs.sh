#!/bin/bash
# =============================================================================
# 00_prepare_inputs.sh
#
# Downloads raw reads from SRA, builds the nf-core/atacseq sample sheet and
# prepares the reference files (genome FASTA, GTF, BWA index, gene BED,
# blacklist).
#
# Manual inputs (download before running):
#   - SraRunTable.csv           SRA Run Selector metadata for the BioProject
#   - <GSE>_series_matrix.txt   GEO series matrix (maps GSM -> sample name)
#
# Sample names are expected to follow <group>_REP<n> (e.g. WT_Skin_REP1).
#
# Requirements: sra-tools (prefetch, fasterq-dump), GNU awk (gawk), gzip,
#               wget, bwa, samtools
# =============================================================================

set -euo pipefail

## ---- editable parameters ---------------------------------------------------
BASE=/path/to/project
DOWNLOAD_DIR=${BASE}/download
FASTQ_DIR=${DOWNLOAD_DIR}/fastq
REF=${BASE}/references

SRA_RUN_TABLE=${DOWNLOAD_DIR}/SraRunTable.csv
SERIES_MATRIX=${DOWNLOAD_DIR}/GSE277902_series_matrix.txt
SAMPLE_FILTER="Skin"     # regex on the GEO sample name; "" keeps all samples
                         # (all must match <group>_REP<n>; suffixes such as _REP1_a are not handled)
SAMPLESHEET=${DOWNLOAD_DIR}/samplesheet.csv

THREADS=16

ENSEMBL_RELEASE=95
ENSEMBL_URL=https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}
FASTA_NAME=Mus_musculus.GRCm38.dna_sm.toplevel.fa
GTF_NAME=Mus_musculus.GRCm38.${ENSEMBL_RELEASE}.gtf
BLACKLIST_URL=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
BLACKLIST_NAME=mm10-blacklist.GRCm38.bed
## ----------------------------------------------------------------------------

mkdir -p "${DOWNLOAD_DIR}/sra_files" "${FASTQ_DIR}" "${REF}/bwa"
cd "${DOWNLOAD_DIR}"

## ---- 1. Run -> GSM -> sample name ------------------------------------------
# Columns are located by header name, as their position varies between exports.
LC_ALL=C gawk -v FPAT='([^,]*)|("[^"]*")' -v OFS='\t' '
NR == 1 {
  sub(/^\xEF\xBB\xBF/, "", $1)
  for (i = 1; i <= NF; i++) { gsub(/"/, "", $i); col[$i] = i }
  if (!("Run" in col) || !("Sample Name" in col)) {
    print "SraRunTable.csv lacks a Run or Sample Name column" > "/dev/stderr"; exit 1
  }
  next
}
{ r = $col["Run"]; s = $col["Sample Name"]; gsub(/"/, "", r); gsub(/"/, "", s); print r, s }
' "${SRA_RUN_TABLE}" > srr_to_gsm.tsv

gawk -F'\t' -v OFS='\t' '
$1 ~ /^!Sample_geo_accession/ { for (i = 2; i <= NF; i++) { gsub(/"/, "", $i); gsm[i] = $i } }
$1 ~ /^!Sample_title/        { for (i = 2; i <= NF; i++) { gsub(/"/, "", $i); name[i] = $i } }
END { for (i in gsm) print gsm[i], name[i] }
' "${SERIES_MATRIX}" > gsm_to_name.tsv

gawk -F'\t' -v OFS='\t' -v filter="${SAMPLE_FILTER}" '
FNR == NR { name[$1] = $2; next }
($2 in name) && name[$2] ~ filter { print $1, $2, name[$2] }
' gsm_to_name.tsv srr_to_gsm.tsv | sort -k3,3 -k1,1 > srr_to_sample.tsv

if [ ! -s srr_to_sample.tsv ]; then
  echo "No runs matched SAMPLE_FILTER='${SAMPLE_FILTER}'" >&2; exit 1
fi
cut -f1 srr_to_sample.tsv > srr_list.txt

## ---- 2. Download and convert to gzipped FASTQ -----------------------------
prefetch --output-directory sra_files --option-file srr_list.txt

# On a cluster, submit one job per run instead.
while read -r srr; do
  [ -s "${FASTQ_DIR}/${srr}_1.fastq.gz" ] && continue   # already converted
  fasterq-dump "sra_files/${srr}/${srr}.sra" \
    --split-files --threads "${THREADS}" --outdir "${FASTQ_DIR}"
  gzip "${FASTQ_DIR}/${srr}"_*.fastq
done < srr_list.txt

## ---- 3. nf-core/atacseq v2.0 sample sheet -----------------------------------
# sample = experimental group, replicate = biological replicate. Runs from the
# same GSM share both values and are merged as technical replicates.
gawk -F'\t' -v dir="${FASTQ_DIR}" '
BEGIN { print "sample,fastq_1,fastq_2,replicate" }
{
  srr = $1; name = $3
  if (!match(name, /^(.+)_REP([0-9]+)$/, m)) {
    print "Cannot parse <group>_REP<n> from sample name: " name > "/dev/stderr"; exit 1
  }
  print m[1] "," dir "/" srr "_1.fastq.gz," dir "/" srr "_2.fastq.gz," m[2]
}' srr_to_sample.tsv > "${SAMPLESHEET}"

## ---- 4. References (Ensembl GRCm38, matching release for FASTA and GTF) -----
cd "${REF}"

[ -s "${FASTA_NAME}" ] || {
  wget -nc "${ENSEMBL_URL}/fasta/mus_musculus/dna/${FASTA_NAME}.gz"
  gunzip "${FASTA_NAME}.gz"
}
[ -s "${GTF_NAME}" ] || {
  wget -nc "${ENSEMBL_URL}/gtf/mus_musculus/${GTF_NAME}.gz"
  gunzip "${GTF_NAME}.gz"
}

# FASTA index, required by the GC-content step in script 03
[ -s "${FASTA_NAME}.fai" ] || samtools faidx "${FASTA_NAME}"

# BWA index, written to references/bwa/ (passed to nf-core as --bwa_index)
[ -s "bwa/${FASTA_NAME}.bwt" ] || bwa index -a bwtsw -p "bwa/${FASTA_NAME}" "${FASTA_NAME}"

# Gene BED6 (0-based start, gene_id, strand) from the GTF
gawk -F'\t' -v OFS='\t' '$3 == "gene" {
  match($9, /gene_id "([^"]+)"/, g)
  print $1, $4 - 1, $5, g[1], 0, $7
}' "${GTF_NAME}" > "${GTF_NAME%.gtf}.bed"

# ENCODE blacklist v2 (mm10), converted from UCSC to Ensembl chromosome names
wget -nc "${BLACKLIST_URL}"
gunzip -c "$(basename "${BLACKLIST_URL}")" \
  | sed -e 's/^chrM\([[:space:]]\)/MT\1/' -e 's/^chr//' > "${BLACKLIST_NAME}"
