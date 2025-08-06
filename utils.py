import time
from Bio import Entrez

Entrez.email = "your_email@example.com"

def load_species_list(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        return set(line.strip() for line in f if line.strip())

def get_species_name(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        for line in record.split("\n"):
            if line.startswith("  ORGANISM"):
                return line.replace("  ORGANISM", "").strip()
    except:
        return "Unknown"

def fetch_fasta_sequences(accession_list):
    batch_size = 10
    sequences = []
    for i in range(0, len(accession_list), batch_size):
        ids = ",".join(accession_list[i:i+batch_size])
        try:
            handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
            sequences.append(handle.read())
            time.sleep(3)
        except Exception as e:
            print(f"Error fetching {ids}: {e}")
    return "".join(sequences)

def extract_target_match(fasta_file):
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                return line.split()[0].replace(">", "")
