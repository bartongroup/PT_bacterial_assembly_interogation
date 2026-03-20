#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V

set -euo pipefail

readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation"
readonly PROTEIN_DIR="${PROJECT_DIR}/protein_fastas"
readonly METADATA_DIR="${PROJECT_DIR}/metadata"
readonly QC_TSV="${METADATA_DIR}/protein_faa_qc.tsv"

printf "isolate_id\tprotein_count\thypothetical_count\tpercent_hypothetical\n" \
    > "${QC_TSV}"

shopt -s nullglob

faa_files=("${PROTEIN_DIR}"/*.faa)

if [[ "${#faa_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No .faa files found in ${PROTEIN_DIR}" >&2
    exit 1
fi

for faa in "${faa_files[@]}"; do
    isolate_id="$(basename "${faa}" .faa)"
    protein_count="$(grep -c '^>' "${faa}" || true)"
    hypothetical_count="$(grep '^>' "${faa}" | grep -ic 'hypothetical protein' || true)"

    if [[ "${protein_count}" -gt 0 ]]; then
        percent_hypothetical="$(
            awk -v hyp="${hypothetical_count}" -v total="${protein_count}" \
                'BEGIN { printf "%.2f", (hyp / total) * 100 }'
        )"
    else
        percent_hypothetical="0.00"
    fi

    printf "%s\t%s\t%s\t%s\n" \
        "${isolate_id}" \
        "${protein_count}" \
        "${hypothetical_count}" \
        "${percent_hypothetical}" \
        >> "${QC_TSV}"
done

echo "[INFO] Protein FASTA QC written to: ${QC_TSV}"
