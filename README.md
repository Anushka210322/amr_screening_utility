# AMR-Screener: In-Silico Antimicrobial Resistance Detection

A lightweight, zero-dependency Python command-line utility designed to screen bacterial *gyrA* (Quinolone Resistance-Determining Region) sequences for fluoroquinolone resistance mutations.

---

## Overview

Screening Sanger consensus reads or isolate assemblies for Point Mutations manually in alignment viewers is slow and prone to human error. This tool automates the detection of critical single-nucleotide polymorphisms (SNPs) in the *gyrA* gene of *Escherichia coli*, focusing on clinical resistance hotspots at **Codon 83** (Serine) and **Codon 87** (Aspartate).

## Key Features

- **Zero External Dependencies:** Built entirely with Python's standard library (`argparse`, `csv`). Runs on any machine without installing third-party packages.
- **Automated Codon Translation:** Translates nucleotide triplets into amino acid positions using standard genetic code tables.
- **Quality Control Flagging:** Automatically detects truncated or incomplete reads and flags them as inconclusive.
- **Standardized Reporting:** Exports machine-ready CSV summaries directly usable in clinical surveillance workflows.

---

## Project Structure

```text
amr_screening_utility/
│
├── data/
│   └── sample_isolates.fasta       # Test clinical FASTA sequences
├── src/
│   └── screener.py                 # Core screening and translation pipeline
├── screening_report.csv            # Generated report output
├── .gitignore                      # Excluded cache and temp files
└── README.md                       # Project documentation