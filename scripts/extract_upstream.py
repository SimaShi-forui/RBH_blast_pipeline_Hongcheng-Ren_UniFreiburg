import os
import sys
from Bio import SeqIO
from BCBio import GFF

# ==== Config ====
BASE_DIR = " "
GENOME_FASTA = os.path.join(BASE_DIR, "haloferax_volcanii.fasta")
GFF_FILE = os.path.join(BASE_DIR, "haloferax_volcanii.gff")
GENE_LIST_FILE = os.path.join(BASE_DIR, "gene_list.txt")
OUTPUT_FASTA = os.path.join(BASE_DIR, "upstream_100bp.fasta")
UPSTREAM_LENGTH = 100
DEBUG = False


# ==== Functions ====

def load_genome(fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Genome FASTA file not found: {fasta_path}")
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def parse_gff_for_genes(gff_path, prefix="HVO_"):
    if not os.path.exists(gff_path):
        raise FileNotFoundError(f"GFF file not found: {gff_path}")
    
    gene_coords = {}
    with open(gff_path, "r") as handle:
        for record in GFF.parse(handle):
            for feature in record.features:
                if feature.type == "gene":
                    gene_id = feature.qualifiers.get("locus_tag", [""])[0]
                    if gene_id.startswith(prefix):
                        gene_coords[gene_id] = {
                            "chrom": record.id,
                            "start": int(feature.location.start),
                            "strand": feature.strand
                        }
    return gene_coords


def read_gene_list(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene list file not found: {path}")
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def extract_upstream_seq(genome, coord, length=100):
    chrom = coord["chrom"]
    start = coord["start"]
    strand = coord["strand"]

    if strand == 1:
        upstream_start = max(0, start - length)
        return genome[chrom].seq[upstream_start:start]
    elif strand == -1:
        return genome[chrom].seq[start:start + length].reverse_complement()
    else:
        return None


def write_fasta(sequences, output_path):
    with open(output_path, "w") as f:
        for gene_id, seq in sequences:
            f.write(f">{gene_id}\n{seq}\n")


# ==== Main ====

def main():
    try:
        genome = load_genome(GENOME_FASTA)
        print(f"[✓] Loaded genome with {len(genome)} contigs.")

        gene_coords = parse_gff_for_genes(GFF_FILE)
        print(f"[✓] Parsed {len(gene_coords)} genes from GFF.")

        gene_list = read_gene_list(GENE_LIST_FILE)
        print(f"[✓] Loaded {len(gene_list)} target genes.")

        missing = [g for g in gene_list if g not in gene_coords]
        if missing:
            print(f"[!] Warning: {len(missing)} genes not found in GFF.")
            if DEBUG:
                print("Missing:", ", ".join(missing[:10]))

        extracted = []
        for gene_id in gene_list:
            if gene_id in gene_coords:
                seq = extract_upstream_seq(genome, gene_coords[gene_id], UPSTREAM_LENGTH)
                if seq:
                    extracted.append((gene_id, seq))

        write_fasta(extracted, OUTPUT_FASTA)
        print(f"[✓] Extracted {len(extracted)} upstream regions to: {OUTPUT_FASTA}")

    except Exception as e:
        print(f"[✗] Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
