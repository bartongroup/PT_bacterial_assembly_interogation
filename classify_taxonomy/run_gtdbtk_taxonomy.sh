#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 16
#$ -jc short
#$ -mods l_hard mfree 62G
#$ -adds l_hard h_vmem 62G

set -euo pipefail

readonly INPUT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/contigs"
readonly OUTPUT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/gtdbtk_taxonomy"
readonly GENOME_DIR="${OUTPUT_DIR}/genomes_for_gtdbtk"
readonly GTDBTK_OUT_DIR="${OUTPUT_DIR}/gtdbtk_results"
readonly LOG_DIR="${OUTPUT_DIR}/logs"

readonly GTDBTK_BIN="/cluster/gjb_lab/pthorpe001/conda/envs/gtdbtk-2.3.2/bin/gtdbtk"
readonly GTDBTK_DB="/cluster/gjb_lab/pthorpe001/conda/envs/gtdbtk-2.3.2/share/gtdbtk-2.3.2/db"
readonly THREADS="16"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${GENOME_DIR}"
mkdir -p "${GTDBTK_OUT_DIR}"
mkdir -p "${LOG_DIR}"

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "[ERROR] Input directory not found: ${INPUT_DIR}" >&2
    exit 1
fi

if [[ ! -x "${GTDBTK_BIN}" ]]; then
    echo "[ERROR] GTDB-Tk binary not found or not executable: ${GTDBTK_BIN}" >&2
    exit 1
fi

if [[ ! -d "${GTDBTK_DB}" ]]; then
    echo "[ERROR] GTDB-Tk database directory not found: ${GTDBTK_DB}" >&2
    exit 1
fi

export GTDBTK_DATA_PATH="${GTDBTK_DB}"
export PATH="$(dirname "${GTDBTK_BIN}"):${PATH}"

echo "[INFO] Input directory: ${INPUT_DIR}"
echo "[INFO] Output directory: ${OUTPUT_DIR}"
echo "[INFO] Genome staging directory: ${GENOME_DIR}"
echo "[INFO] GTDB-Tk output directory: ${GTDBTK_OUT_DIR}"
echo "[INFO] GTDB-Tk binary: ${GTDBTK_BIN}"
echo "[INFO] GTDB-Tk database: ${GTDBTK_DB}"
echo "[INFO] Threads: ${THREADS}"

shopt -s nullglob

fasta_files=("${INPUT_DIR}"/*.fasta)

if [[ "${#fasta_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .fasta files found in ${INPUT_DIR}" >&2
    exit 1
fi

echo "[INFO] Found ${#fasta_files[@]} FASTA assemblies"

rm -f "${GENOME_DIR}"/*.fasta

for fasta in "${fasta_files[@]}"; do
    sample_name="$(basename "${fasta}")"
    ln -s "${fasta}" "${GENOME_DIR}/${sample_name}"
done

printf "sample\tinput_fasta\tstaged_fasta\n" > "${OUTPUT_DIR}/gtdbtk_input_manifest.tsv"

for fasta in "${fasta_files[@]}"; do
    sample_id="$(basename "${fasta}" .fasta)"
    printf "%s\t%s\t%s\n" \
        "${sample_id}" \
        "${fasta}" \
        "${GENOME_DIR}/$(basename "${fasta}")" \
        >> "${OUTPUT_DIR}/gtdbtk_input_manifest.tsv"
done

echo "[INFO] Starting GTDB-Tk classify_wf"

"${GTDBTK_BIN}" classify_wf \
    --genome_dir "${GENOME_DIR}" \
    --out_dir "${GTDBTK_OUT_DIR}" \
    --extension "fasta" \
    --cpus "${THREADS}" \
    --mash_db "${GTDBTK_DB}" \
    --keep_intermediates \
    > "${LOG_DIR}/gtdbtk_classify.stdout.log" \
    2> "${LOG_DIR}/gtdbtk_classify.stderr.log"

echo "[INFO] GTDB-Tk finished"

if [[ -f "${GTDBTK_OUT_DIR}/gtdbtk.bac120.summary.tsv" ]]; then
    echo "[INFO] Bacterial summary found: ${GTDBTK_OUT_DIR}/gtdbtk.bac120.summary.tsv"
fi

if [[ -f "${GTDBTK_OUT_DIR}/gtdbtk.ar53.summary.tsv" ]]; then
    echo "[INFO] Archaeal summary found: ${GTDBTK_OUT_DIR}/gtdbtk.ar53.summary.tsv"
fi

echo "[INFO] Main output folder: ${GTDBTK_OUT_DIR}"
echo "[INFO] Input manifest: ${OUTPUT_DIR}/gtdbtk_input_manifest.tsv"
echo "[INFO] Stdout log: ${LOG_DIR}/gtdbtk_classify.stdout.log"
echo "[INFO] Stderr log: ${LOG_DIR}/gtdbtk_classify.stderr.log"
