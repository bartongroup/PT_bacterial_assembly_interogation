#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 16
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G

set -euo pipefail

readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation"
readonly PROTEIN_DIR="${PROJECT_DIR}/protein_fastas"
readonly EGGNOG_OUT_DIR="${PROJECT_DIR}/eggnog_results"
readonly LOG_DIR="${PROJECT_DIR}/logs"
readonly EGGNOG_DB_DIR="/cluster/gjb_lab/pthorpe001/databases/eggnog"
readonly THREADS="16"

mkdir -p "${EGGNOG_OUT_DIR}"
mkdir -p "${LOG_DIR}"

if ! command -v emapper.py >/dev/null 2>&1; then
    echo "[ERROR] emapper.py not found in PATH. Activate eggnog_env first." >&2
    exit 1
fi

if [[ ! -d "${EGGNOG_DB_DIR}" ]]; then
    echo "[ERROR] eggNOG database directory not found: ${EGGNOG_DB_DIR}" >&2
    exit 1
fi

shopt -s nullglob

faa_files=("${PROTEIN_DIR}"/*.faa)

if [[ "${#faa_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .faa files found in ${PROTEIN_DIR}" >&2
    exit 1
fi

for faa in "${faa_files[@]}"; do
    isolate_id="$(basename "${faa}" .faa)"
    isolate_out_dir="${EGGNOG_OUT_DIR}/${isolate_id}"
    mkdir -p "${isolate_out_dir}"

    echo "[INFO] Running eggNOG-mapper for ${isolate_id}"

    emapper.py \
        --input "${faa}" \
        --itype proteins \
        --output "${isolate_id}" \
        --output_dir "${isolate_out_dir}" \
        --data_dir "${EGGNOG_DB_DIR}" \
        --cpu "${THREADS}" \
        --override \
        --excel \
        > "${LOG_DIR}/${isolate_id}.eggnog.stdout.log" \
        2> "${LOG_DIR}/${isolate_id}.eggnog.stderr.log"
done

echo "[INFO] eggNOG-mapper runs submitted serially"
echo "[INFO] Results directory: ${EGGNOG_OUT_DIR}"
