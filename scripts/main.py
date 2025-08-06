import os
import subprocess
import pandas as pd
from configs import *
from utils import *

def run_blast_for_query(query_file):
    query_name = os.path.splitext(os.path.basename(query_file))[0]
    result_dir = os.path.join(RESULT_ROOT, query_name)
    os.makedirs(result_dir, exist_ok=True)

    target_match = extract_target_match(query_file)
    species_list = load_species_list(SPECIES_FILE)

    blast_output = os.path.join(result_dir, "blast_results.txt")
    reverse_query_file = os.path.join(result_dir, "reverse_query.fasta")
    reverse_output = os.path.join(result_dir, "reverse_blast.txt")
    final_output = os.path.join(result_dir, "final_match_check.txt")

    cmd = [
        "blastp",
        "-query", query_file,
        "-db", ARCHAEA_DB,
        "-evalue", "1e-20",
        "-max_target_seqs", "2000",
        "-outfmt", "6 qseqid sseqid evalue bitscore",
        "-out", blast_output
    ]
    subprocess.run(cmd, check=True)

    df = pd.read_csv(blast_output, sep="\t", names=["query", "accession", "evalue", "bitscore"])
    df["species"] = df["accession"].apply(get_species_name)
    df = df[df["species"].isin(species_list)]

    fasta_content = fetch_fasta_sequences(df["accession"].tolist())
    with open(reverse_query_file, "w") as f:
        f.write(fasta_content)

    cmd_rev = [
        "blastp",
        "-query", reverse_query_file,
        "-db", HALOFERAX_DB,
        "-evalue", "1e-5",
        "-outfmt", "6 qseqid sseqid evalue bitscore",
        "-out", reverse_output
    ]
    subprocess.run(cmd_rev, check=True)

    rev_df = pd.read_csv(reverse_output, sep="\t", names=["accession", "match", "evalue", "bitscore"])
    rev_df = rev_df.sort_values(by=["accession", "evalue"]).groupby("accession").first().reset_index()
    rev_df["Match_Correct"] = rev_df["match"].apply(lambda x: 1 if x == target_match else 0)

    merged = rev_df.merge(df, on="accession", how="left")
    merged = merged[["accession", "Match_Correct", "evalue_x"]]
    merged.columns = ["accession", "Match_Correct", "evalue"]
    merged["species"] = merged["accession"].apply(get_species_name)
    merged.to_csv(final_output, sep="\t", index=False)

def main():
    for file in os.listdir(QUERY_DIR):
        if file.endswith(".fasta"):
            full_path = os.path.join(QUERY_DIR, file)
            run_blast_for_query(full_path)

if __name__ == "__main__":
    main()
