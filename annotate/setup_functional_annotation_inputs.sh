#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -V

set -euo pipefail

readonly RAW_CONTIG_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/contigs"
readonly PROJECT_DIR="/home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation"

readonly INPUT_LINKS_DIR="${PROJECT_DIR}/input_links"
readonly PROTEIN_DIR="${PROJECT_DIR}/protein_fastas"
readonly GENOME_DIR="${PROJECT_DIR}/genome_fastas"
readonly METADATA_DIR="${PROJECT_DIR}/metadata"
readonly LOG_DIR="${PROJECT_DIR}/logs"
readonly EGGNOG_DIR="${PROJECT_DIR}/eggnog_results"
readonly DBCAN_DIR="${PROJECT_DIR}/dbcan_results"
readonly ANTISMASH_DIR="${PROJECT_DIR}/antismash_results"
readonly MATRICES_DIR="${PROJECT_DIR}/matrices"
readonly PGPR_DIR="${PROJECT_DIR}/pgpr_screen"

mkdir -p "${PROJECT_DIR}"
mkdir -p "${INPUT_LINKS_DIR}"
mkdir -p "${PROTEIN_DIR}"
mkdir -p "${GENOME_DIR}"
mkdir -p "${METADATA_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${EGGNOG_DIR}"
mkdir -p "${DBCAN_DIR}"
mkdir -p "${ANTISMASH_DIR}"
mkdir -p "${MATRICES_DIR}"
mkdir -p "${PGPR_DIR}"

readonly MANIFEST_TSV="${METADATA_DIR}/isolate_manifest.tsv"
readonly MISSING_TSV="${METADATA_DIR}/missing_or_incomplete_annotations.tsv"

printf "isolate_id\tgenome_fasta\tannotation_dir\tprotein_faa\tffn\tgff\tgbk\ttsv\n" \
    > "${MANIFEST_TSV}"

printf "isolate_id\tissue\tmissing_item\n" > "${MISSING_TSV}"

shopt -s nullglob

fasta_files=("${RAW_CONTIG_DIR}"/*.fasta)

if [[ "${#fasta_files[@]}" -eq 0 ]]; then
    echo "[ERROR] No genome FASTA files found in ${RAW_CONTIG_DIR}" >&2
    exit 1
fi

for genome_fasta in "${fasta_files[@]}"; do
    isolate_id="$(basename "${genome_fasta}" .fasta)"
    annotation_dir="${RAW_CONTIG_DIR}/${isolate_id}_annotation"

    faa_file="${annotation_dir}/${isolate_id}.faa"
    ffn_file="${annotation_dir}/${isolate_id}.ffn"
    gff_file="${annotation_dir}/${isolate_id}.gff"
    gbk_file="${annotation_dir}/${isolate_id}.gbk"
    tsv_file="${annotation_dir}/${isolate_id}.tsv"

    if [[ ! -d "${annotation_dir}" ]]; then
        printf "%s\tannotation_dir\tmissing\n" "${isolate_id}" >> "${MISSING_TSV}"
        continue
    fi

    for required_file in \
        "${faa_file}" \
        "${ffn_file}" \
        "${gff_file}" \
        "${gbk_file}" \
        "${tsv_file}"; do
        if [[ ! -f "${required_file}" ]]; then
            printf "%s\t%s\tmissing\n" \
                "${isolate_id}" \
                "$(basename "${required_file}")" \
                >> "${MISSING_TSV}"
        fi
    done

    if [[ -f "${genome_fasta}" ]]; then
        ln -sfn "${genome_fasta}" "${GENOME_DIR}/${isolate_id}.fasta"
    fi

    if [[ -f "${faa_file}" ]]; then
        ln -sfn "${faa_file}" "${PROTEIN_DIR}/${isolate_id}.faa"
    fi

    ln -sfn "${annotation_dir}" "${INPUT_LINKS_DIR}/${isolate_id}_annotation"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${isolate_id}" \
        "${genome_fasta}" \
        "${annotation_dir}" \
        "${faa_file}" \
        "${ffn_file}" \
        "${gff_file}" \
        "${gbk_file}" \
        "${tsv_file}" \
        >> "${MANIFEST_TSV}"
done

echo "[INFO] Setup complete"
echo "[INFO] Project directory: ${PROJECT_DIR}"
echo "[INFO] Manifest: ${MANIFEST_TSV}"
echo "[INFO] Missing/incomplete report: ${MISSING_TSV}"
echo "[INFO] Protein FASTA symlinks: ${PROTEIN_DIR}"
echo "[INFO] Genome FASTA symlinks: ${GENOME_DIR}"
