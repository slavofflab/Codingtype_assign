"""
Script to load genome fa file.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import sys

def load_genome(gene_file):
        # Check file extension and open accordingly
    if gene_file.endswith(".gz"):
        try:
            gf = gzip.open(gene_file, 'rt')  # Use 'rt' for text mode
        except IOError as e:
            print(f"Error opening gzipped file: {e}")
            sys.exit(1)
    else:
        try:
            gf = open(gene_file, 'r')
        except IOError as e:
            print(f"Error opening file: {e}")
            sys.exit(1)
    chrom_seq = {}

    print("Start loading gene file")
    try:
        for record in SeqIO.parse(gf, "fasta"):
            chrom_seq[record.id] = Seq(str(record.seq)).upper()
            print(f"{record.id} loaded")
    except Exception as e:
        print(f"Error parsing file: {e}")
        sys.exit(1)
    finally:
        gf.close()

    return chrom_seq
