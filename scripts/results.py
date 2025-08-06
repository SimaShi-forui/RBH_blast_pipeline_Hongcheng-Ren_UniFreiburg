
import os
import time
import pandas as pd
from Bio import Entrez
from configs import EMAIL, FINAL_RESULT_DIR
from utils import get_species_name

Entrez.email = EMAIL

def get_species_and_locus(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        species, locus = "Unknown", ""
        for line in record.splitlines():
            if line.startswith("  ORGANISM"):
                species = " ".join(line.split()[1:3])
            if "/locus_tag=" in line:
                locus = line.split("=")[1].replace('"', '')
        return species, locus
    except:
        return "Unknown", ""

def merge_match_results():
    files = [f for f in os.listdir(FINAL_RESULT_DIR) if f.endswith(".txt") and "_final_match_check_" in f]
    all_accessions = set()
    dataframes = {}

    for f in files:
        df = pd.read_csv(os.path.join(FINAL_RESULT_DIR, f), sep="\t")
        df = df[["accession", "Match_Correct"]]
        dataframes[f] = df
        all_accessions.update(df["accession"])

    acc_info = {}
    for acc in all_accessions:
        acc_info[acc] = get_species_and_locus(acc)
        time.sleep(0.4)

    records = []
    for acc in all_accessions:
        record = {
            "Accession": acc,
            "Species": acc_info[acc][0],
            "Locus": acc_info[acc][1]
        }
        for f in files:
            match = dataframes[f][dataframes[f]["accession"] == acc]
            record[f] = int(match["Match_Correct"].values[0]) if not match.empty else -1
        records.append(record)

    merged_df = pd.DataFrame(records)
    output_file = os.path.join(FINAL_RESULT_DIR, "merged_match_results.tsv")
    merged_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    merge_match_results()
