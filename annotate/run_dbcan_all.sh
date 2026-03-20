#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 8
#$ -jc long
#$ -mods l_hard mfree 16G
#$ -adds l_hard h_vmem 16G

set -euo pipefail

readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation"
readonly PROTEIN_DIR="${PROJECT_DIR}/protein_fastas"
readonly DBCAN_OUT_DIR="${PROJECT_DIR}/dbcan_results"
readonly METADATA_DIR="${PROJECT_DIR}/metadata"
readonly LOG_DIR="${PROJECT_DIR}/logs"

readonly DBCAN_DB_DIR="/cluster/gjb_lab/pthorpe001/databases/dbcan_db"
readonly THREADS="8"

mkdir -p "${DBCAN_OUT_DIR}"
mkdir -p "${METADATA_DIR}"
mkdir -p "${LOG_DIR}"

if ! command -v run_dbcan >/dev/null 2>&1; then
    echo "[ERROR] run_dbcan not found in PATH. Activate dbcan_env first." >&2
    exit 1
fi

if [[ ! -d "${DBCAN_DB_DIR}" ]]; then
    echo "[ERROR] dbCAN database directory not found: ${DBCAN_DB_DIR}" >&2
    exit 1
fi

shopt -s nullglob

faa_files=("${PROTEIN_DIR}"/*.faa)

if [[ "${#faa_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .faa files found in ${PROTEIN_DIR}" >&2
    exit 1
fi

readonly INPUT_MANIFEST="${METADATA_DIR}/dbcan_input_manifest.tsv"

printf "isolate_id\tprotein_faa\toutput_dir\n" > "${INPUT_MANIFEST}"

for faa in "${faa_files[@]}"; do
    isolate_id="$(basename "${faa}" .faa)"
    isolate_out_dir="${DBCAN_OUT_DIR}/${isolate_id}"

    mkdir -p "${isolate_out_dir}"

    printf "%s\t%s\t%s\n" \
        "${isolate_id}" \
        "${faa}" \
        "${isolate_out_dir}" \
        >> "${INPUT_MANIFEST}"

    if [[ -f "${isolate_out_dir}/overview.txt" ]] || \
       [[ -f "${isolate_out_dir}/dbCAN-sub.hmm.out" ]] || \
       [[ -f "${isolate_out_dir}/diamond.out" ]]; then
        echo "[INFO] dbCAN output already present for ${isolate_id}, skipping"
        continue
    fi

    echo "[INFO] Running dbCAN for ${isolate_id}"

    run_dbcan CAZyme_annotation \
        --input_raw_data "${faa}" \
        --mode protein \
        --output_dir "${isolate_out_dir}" \
        --db_dir "${DBCAN_DB_DIR}" \
        > "${LOG_DIR}/${isolate_id}.dbcan.stdout.log" \
        2> "${LOG_DIR}/${isolate_id}.dbcan.stderr.log"
done

echo "[INFO] dbCAN run complete"
echo "[INFO] Results directory: ${DBCAN_OUT_DIR}"
echo "[INFO] Manifest: ${INPUT_MANIFEST}"
