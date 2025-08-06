# RBH-based BLAST Pipeline for Archaeal Protein Homolog Detection

This repository provides a Python-based pipeline for identifying homologous proteins in archaea using reciprocal BLAST (RBH), and constructing phylogenetic trees based on sequence similarity.

---

**Author**: Hongcheng Ren  
**Lab**: Molecular Biology of Archaea 
**Affiliation**: Institute for Biology II- Microbiology, University of Freiburg, Germany  
**Contact**: simashi.orbit (at) gmail.com  
**Repository Maintainer**: [SimaShi-forui](https://github.com/SimaShi-forui)

---

## Project Features

- Automatic parsing of multiple `.fasta` query files
- Reciprocal BLAST to identify orthologous candidates
- Species filtering via Entrez
- Reverse BLAST confirmation
- Summary matrix for all matched proteins
- Phylogenetic tree construction using FastTree

---

## Repository Structure
.
├── main.py # Batch process all RBH analysis
├── configs.py # Centralized configuration
├── utils.py # Helper functions
├── enrich/
│ ├── merge_match_results.py
│ ├── summarize_final_match.py
│ └── run_fasttree.py
├── results/ # Auto-generated results
├── data/queries/ # Place your query .fasta files here
├── .gitignore
├── README.md
└── LICENSE

## Dependencies
Biopython
pandas
blast+
FastTree
Use the following to install:

pip install biopython pandas

## Quick Start

1. Place your `.fasta` files in `data/queries/`
2. Set your species list in `data/species_list.txt`
3. Run the main script:

```bash
python main.py
python enrich/merge_match_results.py
python enrich/summarize_final_match.py
python enrich/run_fasttree.py

