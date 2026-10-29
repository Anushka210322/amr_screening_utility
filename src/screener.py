"""
AMR Screening Utility
Screens bacterial gyrA sequences for fluoroquinolone resistance mutations.
Pure Python - Uses only standard built-ins (argparse, csv, os, sys).
"""

import argparse
import csv
import os
import sys

# Standard Genetic Code Dictionary
CODON_MAP = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*", "TGA": "*",
    "TGC": "C", "TGT": "C", "TGG": "W",
}

# Known resistance loci in E. coli gyrA (QRDR)
RESISTANCE_RULES = {
    83: {"WT": "S", "resistant_alleles": {"L": "Ser83Leu", "W": "Ser83Trp"}},
    87: {"WT": "D", "resistant_alleles": {"N": "Asp87Asn", "Y": "Asp87Tyr"}},
}


def parse_fasta(filepath):
    """Parses multi-line FASTA files and cleans stray gaps/spaces."""
    header = None
    seq_chunks = []
    clean_table = str.maketrans("", "", "-. \t\r\n")

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).translate(clean_table).upper()
                header = line_str[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line_str)
        if header is not None:
            yield header, "".join(seq_chunks).translate(clean_table).upper()


def extract_codon_and_aa(dna_seq, aa_position):
    """Extracts a 3-base codon at 1-based amino acid position and translates it."""
    start_idx = (aa_position - 1) * 3
    end_idx = start_idx + 3

    if len(dna_seq) < end_idx:
        return "---", "TRUNCATED"

    codon = dna_seq[start_idx:end_idx]
    aa = CODON_MAP.get(codon, "Unknown")
    return codon, aa


def screen_isolate(sample_id, dna_seq):
    """Performs codon extraction, mutation identification, and phenotype assignment."""
    row = {
        "sample_id": sample_id,
        "seq_len_bp": len(dna_seq),
    }

    mutations = []
    is_resistant = False
    qc_failed = False

    for pos, info in RESISTANCE_RULES.items():
        codon, aa = extract_codon_and_aa(dna_seq, pos)
        wt = info["WT"]

        row[f"pos{pos}_codon"] = codon
        row[f"pos{pos}_aa"] = aa

        if aa == "TRUNCATED":
            row[f"pos{pos}_mutation"] = "QC_FAIL_TOO_SHORT"
            qc_failed = True
        elif aa == wt:
            row[f"pos{pos}_mutation"] = "WT"
        else:
            mut_name = info["resistant_alleles"].get(aa, f"{wt}{pos}{aa}")
            row[f"pos{pos}_mutation"] = mut_name
            mutations.append(mut_name)
            is_resistant = True

    if qc_failed:
        row["clinical_phenotype"] = "INCONCLUSIVE (Sequence Truncated)"
    elif is_resistant:
        row["clinical_phenotype"] = f"RESISTANT ({', '.join(mutations)})"
    else:
        row["clinical_phenotype"] = "SENSITIVE (Susceptible)"

    return row


def run(input_path, output_path):
    """Coordinates file reading, processing, and CSV output."""
    if not os.path.exists(input_path):
        print(f"Error: Input file '{input_path}' not found.")
        sys.exit(1)

    records = list(parse_fasta(input_path))
    if not records:
        print("Error: No valid FASTA sequences found in the input file.")
        sys.exit(1)

    screened_data = [screen_isolate(header, seq) for header, seq in records]

    fieldnames = list(screened_data[0].keys())
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(screened_data)

    print("-" * 55)
    print(f" Successfully screened {len(screened_data)} clinical isolates.")
    print(f" Report generated: {output_path}")
    print("-" * 55)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CLI tool to screen bacterial FASTA sequences for AMR mutations."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to input FASTA file"
    )
    parser.add_argument(
        "-o", "--output", default="amr_report.csv", help="Path to output CSV file"
    )

    args = parser.parse_args()
    run(args.input, args.output)