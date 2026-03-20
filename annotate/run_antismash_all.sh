#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 16
#$ -jc long
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G

set -euo pipefail

readonly RAW_CONTIG_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/contigs"
readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation"
readonly ANTISMASH_OUT_DIR="${PROJECT_DIR}/antismash_results"
readonly LOG_DIR="${PROJECT_DIR}/logs"

readonly ANTISMASH_DB="/cluster/gjb_lab/pthorpe001/conda/envs/antismash/lib/python3.10/site-packages/antismash/databases"
readonly THREADS="16"

mkdir -p "${ANTISMASH_OUT_DIR}"
mkdir -p "${LOG_DIR}"

if ! command -v antismash >/dev/null 2>&1; then
    echo "[ERROR] antismash not found in PATH. Activate the antismash environment first." >&2
    exit 1
fi

if [[ ! -d "${ANTISMASH_DB}" ]]; then
    echo "[ERROR] antiSMASH database directory not found: ${ANTISMASH_DB}" >&2
    exit 1
fi

shopt -s nullglob

gbk_files=("${RAW_CONTIG_DIR}"/*_annotation/*.gbk)

if [[ "${#gbk_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .gbk files found under ${RAW_CONTIG_DIR}/*_annotation/" >&2
    exit 1
fi

printf "isolate_id\tgbk_input\toutput_dir\n" \
    > "${PROJECT_DIR}/metadata/antismash_input_manifest.tsv"

for gbk in "${gbk_files[@]}"; do
    isolate_id="$(basename "${gbk}" .gbk)"
    isolate_out_dir="${ANTISMASH_OUT_DIR}/${isolate_id}"

    mkdir -p "${isolate_out_dir}"

    printf "%s\t%s\t%s\n" \
        "${isolate_id}" \
        "${gbk}" \
        "${isolate_out_dir}" \
        >> "${PROJECT_DIR}/metadata/antismash_input_manifest.tsv"

    if [[ -f "${isolate_out_dir}/${isolate_id}.zip" ]] || \
       [[ -f "${isolate_out_dir}/index.html" ]]; then
        echo "[INFO] antiSMASH output already present for ${isolate_id}, skipping"
        continue
    fi

    echo "[INFO] Running antiSMASH for ${isolate_id}"

    antismash \
        --output-dir "${isolate_out_dir}" \
        --output-basename "${isolate_id}" \
        --taxon bacteria \
        --pfam2go \
        --cpus "${THREADS}" \
        --databases "${ANTISMASH_DB}" \
        "${gbk}" \
        > "${LOG_DIR}/${isolate_id}.antismash.stdout.log" \
        2> "${LOG_DIR}/${isolate_id}.antismash.stderr.log"
done

echo "[INFO] antiSMASH run complete"
echo "[INFO] Results directory: ${ANTISMASH_OUT_DIR}"
echo "[INFO] Manifest: ${PROJECT_DIR}/metadata/antismash_input_manifest.tsv"
