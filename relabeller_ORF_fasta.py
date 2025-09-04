#### This script requires biopython
from Bio import SeqIO
import argparse

def rename_fasta_headers(input_file, output_file, prefix='seq'):
    """
    Rename all headers in a FASTA file sequentially.
    
    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
        prefix (str, optional): Prefix for the new headers. Defaults to 'seq'.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    
    for i, record in enumerate(records, 1):
        record.id = f"{prefix}_{i}"
        record.description = record.id
    
    SeqIO.write(records, output_file, "fasta")
    print(f"Renamed {len(records)} sequences. Output written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers sequentially.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-p", "--prefix", default="seq", help="Prefix for new headers (default: seq)")
    
    args = parser.parse_args()
    
    rename_fasta_headers(args.input, args.output, args.prefix)

if __name__ == "__main__":
    main()
