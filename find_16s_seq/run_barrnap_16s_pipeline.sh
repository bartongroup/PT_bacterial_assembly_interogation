#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V

set -euo pipefail



readonly INPUT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/contigs"
readonly OUTPUT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/find_16s_sequences"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly GFF_DIR="${OUTPUT_DIR}/gff"
readonly FASTA_DIR="${OUTPUT_DIR}/16s_fasta"
readonly SUMMARY_DIR="${OUTPUT_DIR}/summaries"
readonly LOG_DIR="${OUTPUT_DIR}/logs"
readonly PYTHON_SCRIPT="${SCRIPT_DIR}/extract_16s_from_barrnap.py"

mkdir -p "${GFF_DIR}" "${FASTA_DIR}" "${SUMMARY_DIR}" "${LOG_DIR}"


if ! command -v barrnap >/dev/null 2>&1; then
    echo "[ERROR] barrnap not found in PATH." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 not found in PATH." >&2
    exit 1
fi

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "[ERROR] Python script not found: ${PYTHON_SCRIPT}" >&2
    exit 1
fi


shopt -s nullglob

fasta_files=("${INPUT_DIR}"/*.fasta)

if [[ "${#fasta_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .fasta files found in ${INPUT_DIR}" >&2
    exit 1
fi

for fasta in "${fasta_files[@]}"; do
    isolate="$(basename "${fasta}" .fasta)"
    gff_out="${GFF_DIR}/${isolate}.barrnap.gff"
    log_out="${LOG_DIR}/${isolate}.barrnap.log"
    fasta_out="${FASTA_DIR}/${isolate}.16s.fasta"
    summary_out="${SUMMARY_DIR}/${isolate}.16s_summary.tsv"

    echo "[INFO] Running barrnap for ${isolate}"

    barrnap \
        --kingdom bac \
        --threads 1 \
        --evalue 1e-06 \
        "${fasta}" \
        > "${gff_out}" \
        2> "${log_out}"

    echo "[INFO] Extracting 16S sequences for ${isolate}"

    python3 "${PYTHON_SCRIPT}" \
        --assembly "${fasta}" \
        --gff "${gff_out}" \
        --isolate "${isolate}" \
        --out-fasta "${fasta_out}" \
        --out-summary "${summary_out}"
done

echo "[INFO] Building combined outputs"

python3 "${PYTHON_SCRIPT}" \
    --combine-fasta-dir "${FASTA_DIR}" \
    --combine-summary-dir "${SUMMARY_DIR}" \
    --combined-fasta "${FASTA_DIR}/all_isolates_16s.fasta" \
    --combined-summary "${SUMMARY_DIR}/all_isolates_16s_summary.tsv"

echo "[INFO] Done"
echo "[INFO] GFF files: ${GFF_DIR}"
echo "[INFO] Per-isolate 16S FASTA files: ${FASTA_DIR}"
echo "[INFO] Per-isolate summary files: ${SUMMARY_DIR}"
echo "[INFO] Combined FASTA: ${FASTA_DIR}/all_isolates_16s.fasta"
echo "[INFO] Combined summary: ${SUMMARY_DIR}/all_isolates_16s_summary.tsv"
