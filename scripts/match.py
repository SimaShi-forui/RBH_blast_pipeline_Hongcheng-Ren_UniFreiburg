
import os
import time
import pandas as pd
from Bio import Entrez
from configs import EMAIL, FINAL_RESULT_DIR, SUMMARY_OUTPUT_DIR
from utils import get_species_name

Entrez.email = EMAIL
os.makedirs(SUMMARY_OUTPUT_DIR, exist_ok=True)

def fetch_info(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        species = "Unknown"
        locus = "Unknown"
        for line in record.splitlines():
            if line.startswith("  ORGANISM"):
                species = " ".join(line.split()[1:3])
            if "/locus_tag=" in line or "/gene=" in line:
                locus = line.split("=")[1].strip('"')
        return species, locus
    except Exception as e:
        return "Unknown", "Unknown"

def summarize_match_tables():
    summary_dict = {}

    for filename in os.listdir(FINAL_RESULT_DIR):
        if filename.endswith(".txt") and "final_match_check" in filename:
            df = pd.read_csv(os.path.join(FINAL_RESULT_DIR, filename), sep="\t")
            species_list, locus_list = [], []
            for acc in df["accession"]:
                species, locus = fetch_info(acc)
                species_list.append(species)
                locus_list.append(locus)
                time.sleep(0.4)
            df["species_full"] = species_list
            df["locus_tag"] = locus_list

            out_file = os.path.join(SUMMARY_OUTPUT_DIR, filename)
            df.to_csv(out_file, sep="\t", index=False)

            protein_tag = filename.replace("_final_match_check_", "_").replace(".txt", "")
            for _, row in df.iterrows():
                sp = row["species_full"]
                val = row["Match_Correct"]
                if sp not in summary_dict:
                    summary_dict[sp] = {}
                summary_dict[sp][protein_tag] = val

    summary_df = pd.DataFrame.from_dict(summary_dict, orient="index").fillna(-1).astype(int)
    summary_df.index.name = "Species"
    summary_df.to_csv(os.path.join(SUMMARY_OUTPUT_DIR, "summary_matrix.tsv"), sep="\t")

if __name__ == "__main__":
    summarize_match_tables()
