import os

RESULT_ROOT = "results"
QUERY_DIR = "data/queries"
SPECIES_FILE = "data/species_list.txt"
ARCHAEA_DB = r"C:\\blast\\refseqselect_db"
HALOFERAX_DB = r"C:\\blast\\haloferax_db"
EMAIL = "your_email@example.com"  # masked for GitHub

FINAL_RESULT_DIR = os.path.join(RESULT_ROOT, "finals")
SUMMARY_OUTPUT_DIR = os.path.join(RESULT_ROOT, "enriched")
FASTTREE_EXE = r"C:\\FastTree\\FastTree.exe"
FASTTREE_INPUT = os.path.join(RESULT_ROOT, "alignment", "with_gap_output.fasta")
FASTTREE_OUTPUT = os.path.join(RESULT_ROOT, "tree", "phylogenetic_tree.nwk")
