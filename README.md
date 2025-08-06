
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


├── data/ # Data directory (input/output files)
│ ├── output/ # BLAST and result files
│ └── tree/ # Phylogenetic tree output files
│
├── scripts/ # All pipeline-related scripts
│ ├── main.py # Main pipeline execution script
│ ├── configs.py # Configuration and parameters
│ ├── utils.py # Helper functions
│ ├── match.py # Reciprocal BLAST matching
│ ├── results.py # Post-processing and result filtering
│ ├── fasttree.py # Tree construction with FastTree
│ └── extract_upstream.py # Upstream region extractor from GFF
│
├── LICENSE # License information (MIT)
├── README.md # Project overview and documentation


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
python results.py
python match.py
python fasttree.py

