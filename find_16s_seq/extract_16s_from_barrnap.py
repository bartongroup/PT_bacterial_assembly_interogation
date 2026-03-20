#!/usr/bin/env python3

"""
Extract 16S rRNA sequences from barrnap GFF output.

This script supports two modes:

1. Per-isolate extraction mode
   Reads an assembly FASTA and its barrnap GFF file, extracts all features
   labelled as 16S rRNA, writes a per-isolate FASTA file, and writes a
   per-isolate summary TSV.

2. Combine mode
   Concatenates all per-isolate 16S FASTA files into one combined FASTA and
   concatenates all per-isolate summary TSV files into one combined summary TSV.

The script is intended for bacterial isolate assemblies.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO
from Bio.Seq import Seq

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract 16S sequences from barrnap GFF output."
    )

    parser.add_argument(
        "--assembly",
        required=False,
        help="Input assembly FASTA file.",
    )
    parser.add_argument(
        "--gff",
        required=False,
        help="Input barrnap GFF file.",
    )
    parser.add_argument(
        "--isolate",
        required=False,
        help="Isolate name.",
    )
    parser.add_argument(
        "--out-fasta",
        required=False,
        help="Output FASTA file containing extracted 16S sequences.",
    )
    parser.add_argument(
        "--out-summary",
        required=False,
        help="Output TSV summary file for extracted 16S sequences.",
    )

    parser.add_argument(
        "--combine-fasta-dir",
        required=False,
        help="Directory containing per-isolate 16S FASTA files to combine.",
    )
    parser.add_argument(
        "--combine-summary-dir",
        required=False,
        help="Directory containing per-isolate 16S summary TSV files to combine.",
    )
    parser.add_argument(
        "--combined-fasta",
        required=False,
        help="Output combined FASTA file.",
    )
    parser.add_argument(
        "--combined-summary",
        required=False,
        help="Output combined summary TSV file.",
    )

    return parser.parse_args()


def load_assembly_sequences(assembly_path: Path) -> Dict[str, str]:
    """
    Load assembly sequences from a FASTA file.

    Parameters
    ----------
    assembly_path : Path
        Path to the input assembly FASTA file.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping contig identifiers to nucleotide sequence strings.
    """
    sequences: Dict[str, str] = {}

    for record in SeqIO.parse(handle=str(assembly_path), format="fasta"):
        sequences[record.id] = str(record.seq)

    return sequences


def is_16s_feature(feature_type: str, attributes: str) -> bool:
    """
    Determine whether a GFF feature corresponds to a 16S rRNA.

    Parameters
    ----------
    feature_type : str
        GFF feature type.
    attributes : str
        GFF attributes field.

    Returns
    -------
    bool
        True if the feature appears to be a 16S rRNA, otherwise False.
    """
    if feature_type != "rRNA":
        return False

    return "16S" in attributes


def extract_16s_records(
    assembly_path: Path,
    gff_path: Path,
    isolate_name: str,
) -> List[dict]:
    """
    Extract 16S records from an assembly using barrnap GFF annotations.

    Parameters
    ----------
    assembly_path : Path
        Path to the assembly FASTA file.
    gff_path : Path
        Path to the barrnap GFF file.
    isolate_name : str
        Name of the isolate.

    Returns
    -------
    List[dict]
        List of extracted 16S records, each containing metadata and sequence.
    """
    sequences = load_assembly_sequences(assembly_path=assembly_path)
    records: List[dict] = []

    with gff_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            contig = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            if not is_16s_feature(
                feature_type=feature_type,
                attributes=attributes,
            ):
                continue

            if contig not in sequences:
                raise KeyError(
                    f"Contig '{contig}' from GFF not found in assembly FASTA "
                    f"for isolate '{isolate_name}'."
                )

            contig_seq = sequences[contig]

            sequence = contig_seq[start - 1:end]

            if strand == "-":
                sequence = str(Seq(sequence).reverse_complement())
            else:
                sequence = str(sequence)

            length = len(sequence)

            record_id = (
                f"{isolate_name}|{contig}:{start}-{end}|{strand}|16S_rRNA"
            )

            records.append(
                {
                    "isolate": isolate_name,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": length,
                    "product": "16S_rRNA",
                    "attributes": attributes,
                    "record_id": record_id,
                    "sequence": sequence,
                }
            )

    return records


def write_fasta(records: List[dict], output_path: Path) -> None:
    """
    Write extracted 16S records to a FASTA file.

    Parameters
    ----------
    records : List[dict]
        Extracted 16S records.
    output_path : Path
        Output FASTA file path.
    """
    with output_path.open(mode="w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record['record_id']}\n")
            sequence = str(record["sequence"])
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")


def write_summary(records: List[dict], output_path: Path) -> None:
    """
    Write extracted 16S records to a TSV summary file.

    Parameters
    ----------
    records : List[dict]
        Extracted 16S records.
    output_path : Path
        Output TSV file path.
    """
    header = [
        "isolate",
        "contig",
        "start",
        "end",
        "strand",
        "length",
        "product",
        "attributes",
    ]

    with output_path.open(mode="w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")

        for record in records:
            row = [
                str(record["isolate"]),
                str(record["contig"]),
                str(record["start"]),
                str(record["end"]),
                str(record["strand"]),
                str(record["length"]),
                str(record["product"]),
                str(record["attributes"]),
            ]
            handle.write("\t".join(row) + "\n")


def combine_fastas(input_dir: Path, output_path: Path) -> None:
    """
    Combine per-isolate FASTA files into one FASTA file.

    Parameters
    ----------
    input_dir : Path
        Directory containing per-isolate FASTA files.
    output_path : Path
        Output combined FASTA file.
    """
    fasta_files = sorted(input_dir.glob("*.16s.fasta"))

    with output_path.open(mode="w", encoding="utf-8") as out_handle:
        for fasta_file in fasta_files:
            if fasta_file.resolve() == output_path.resolve():
                continue

            with fasta_file.open(mode="r", encoding="utf-8") as in_handle:
                for line in in_handle:
                    out_handle.write(line)


def combine_summaries(input_dir: Path, output_path: Path) -> None:
    """
    Combine per-isolate summary TSV files into one summary TSV file.

    Parameters
    ----------
    input_dir : Path
        Directory containing per-isolate summary TSV files.
    output_path : Path
        Output combined summary TSV file.
    """
    summary_files = sorted(input_dir.glob("*.16s_summary.tsv"))
    header_written = False

    with output_path.open(mode="w", encoding="utf-8") as out_handle:
        for summary_file in summary_files:
            if summary_file.resolve() == output_path.resolve():
                continue

            with summary_file.open(mode="r", encoding="utf-8") as in_handle:
                for line_number, line in enumerate(in_handle):
                    if line_number == 0:
                        if not header_written:
                            out_handle.write(line)
                            header_written = True
                    else:
                        out_handle.write(line)


def run_extract_mode(args: argparse.Namespace) -> None:
    """
    Run per-isolate extraction mode.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    required = [
        args.assembly,
        args.gff,
        args.isolate,
        args.out_fasta,
        args.out_summary,
    ]

    if any(value is None for value in required):
        raise ValueError(
            "Extraction mode requires --assembly, --gff, --isolate, "
            "--out-fasta, and --out-summary."
        )

    records = extract_16s_records(
        assembly_path=Path(args.assembly),
        gff_path=Path(args.gff),
        isolate_name=args.isolate,
    )

    write_fasta(records=records, output_path=Path(args.out_fasta))
    write_summary(records=records, output_path=Path(args.out_summary))


def run_combine_mode(args: argparse.Namespace) -> None:
    """
    Run combined output mode.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    required = [
        args.combine_fasta_dir,
        args.combine_summary_dir,
        args.combined_fasta,
        args.combined_summary,
    ]

    if any(value is None for value in required):
        raise ValueError(
            "Combine mode requires --combine-fasta-dir, "
            "--combine-summary-dir, --combined-fasta, and "
            "--combined-summary."
        )

    combine_fastas(
        input_dir=Path(args.combine_fasta_dir),
        output_path=Path(args.combined_fasta),
    )
    combine_summaries(
        input_dir=Path(args.combine_summary_dir),
        output_path=Path(args.combined_summary),
    )


def main() -> None:
    """
    Run the script.
    """
    args = parse_args()

    extract_mode_requested = any(
        value is not None
        for value in [
            args.assembly,
            args.gff,
            args.isolate,
            args.out_fasta,
            args.out_summary,
        ]
    )

    combine_mode_requested = any(
        value is not None
        for value in [
            args.combine_fasta_dir,
            args.combine_summary_dir,
            args.combined_fasta,
            args.combined_summary,
        ]
    )

    if extract_mode_requested and combine_mode_requested:
        raise ValueError(
            "Please run either extraction mode or combine mode, not both."
        )

    if extract_mode_requested:
        run_extract_mode(args=args)
        return

    if combine_mode_requested:
        run_combine_mode(args=args)
        return

    raise ValueError(
        "No valid mode selected. Use extraction mode arguments or combine "
        "mode arguments."
    )


if __name__ == "__main__":
    main()