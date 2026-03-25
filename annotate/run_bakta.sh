#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 16
#$ -jc short
#$ -mods l_hard mfree 16G
#$ -adds l_hard h_vmem 16G

# conda activate hifi_assembly

set -euo pipefail

readonly INPUT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation/protein_fastas"
readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/harmonised_annotation"
readonly RESULTS_DIR="${PROJECT_DIR}/bakta_proteins_results"
readonly METADATA_DIR="${PROJECT_DIR}/metadata"
readonly LOG_DIR="${PROJECT_DIR}/logs"
readonly THREADS="16"

readonly BAKTA_DB_DEFAULT="/cluster/gjb_lab/pthorpe001/databases/bakta_db"

mkdir -p "${RESULTS_DIR}" "${METADATA_DIR}" "${LOG_DIR}"

if ! command -v bakta_proteins >/dev/null 2>&1; then
    echo "[ERROR] bakta_proteins not found in PATH. Activate the Bakta environment first." >&2
    exit 1
fi

if [[ -n "${BAKTA_DB:-}" ]]; then
    readonly BAKTA_DB_PATH="${BAKTA_DB}"
else
    readonly BAKTA_DB_PATH="${BAKTA_DB_DEFAULT}"
fi

if [[ ! -d "${BAKTA_DB_PATH}" ]]; then
    echo "[ERROR] Bakta database directory not found: ${BAKTA_DB_PATH}" >&2
    echo "[ERROR] Please export BAKTA_DB or edit BAKTA_DB_DEFAULT in the script." >&2
    exit 1
fi

printf "isolate_id\tprotein_faa\toutput_dir\n" > "${METADATA_DIR}/bakta_proteins_input_manifest.tsv"

shopt -s nullglob
faa_files=("${INPUT_DIR}"/*.faa)

if [[ "${#faa_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .faa files found in ${INPUT_DIR}" >&2
    exit 1
fi

for faa in "${faa_files[@]}"; do
    isolate_id="$(basename "${faa}" .faa)"
    isolate_out_dir="${RESULTS_DIR}/${isolate_id}"

    mkdir -p "${isolate_out_dir}"

    printf "%s\t%s\t%s\n" \
        "${isolate_id}" \
        "${faa}" \
        "${isolate_out_dir}" \
        >> "${METADATA_DIR}/bakta_proteins_input_manifest.tsv"

    if [[ -f "${isolate_out_dir}/${isolate_id}.tsv" ]] || \
       [[ -f "${isolate_out_dir}/${isolate_id}.json" ]]; then
        echo "[INFO] Bakta protein output already present for ${isolate_id}, skipping"
        continue
    fi

    echo "[INFO] Running bakta_proteins for ${isolate_id}"

    bakta_proteins \
        --db "${BAKTA_DB_PATH}" \
        --output "${isolate_out_dir}" \
        --prefix "${isolate_id}" \
        --threads "${THREADS}" \
        --force \
        --db /home/pthorpe001/data/databases/bakta \
        "${faa}" \
        > "${LOG_DIR}/${isolate_id}.bakta_proteins.stdout.log" \
        2> "${LOG_DIR}/${isolate_id}.bakta_proteins.stderr.log"
done

echo "[INFO] bakta_proteins run complete"
echo "[INFO] Results directory: ${RESULTS_DIR}"
echo "[INFO] Manifest: ${METADATA_DIR}/bakta_proteins_input_manifest.tsv"