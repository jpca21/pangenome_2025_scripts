#!/usr/bin/env python3
# LAST VERSION: July 9, 2024 (Converted from Perl to Python3)

import sys
import os
import subprocess
import re
from Bio import SeqIO

# Set default inputs
OG_list = sys.argv[1] if len(sys.argv) > 1 else "Orthogroups_SingleCopyOrthologues.txt"
extensive_list = sys.argv[2] if len(sys.argv) > 2 else "Orthogroups.txt"
extensive_fasta = sys.argv[3] if len(sys.argv) > 3 else "ALL_FAA"
header_list = sys.argv[4] if len(sys.argv) > 4 else "genome_to_id.tsv"
OUTPUT = sys.argv[5] if len(sys.argv) > 5 else "OUTPUT_CONC.ALN.fasta"

# Configuring MAFFT
mafft_exe = "mafft"
options = "--maxiterate 1000 --localpair --op 5"

# Obtaining OG list
print("### Obtaining OG list")
content_OG = []
try:
    with open(OG_list, 'r') as f:
        content_OG = [line.strip() for line in f.readlines()]
except FileNotFoundError:
    print(f"Error: Could not open file {OG_list}")
    sys.exit(1)

# Obtaining Extensive data (as a dictionary)
print("### Obtaining Orthogroups to IDs data")
content_data = {}
try:
    with open(extensive_list, 'r') as f:
        for line in f:
            line = line.strip().replace('\r', '')
            # Using regex to match the pattern
            match = re.match(r'^(\S+):\s*(.+)$', line)
            if match:
                a = match.group(1)
                b = match.group(2)
                content_data[a] = b
except FileNotFoundError:
    print(f"Error: Could not open file {extensive_list}")
    sys.exit(1)

# Obtaining Headers list
print("### Obtaining Genome to IDs data")
content_headers = {}
try:
    with open(header_list, 'r') as f:
        for line in f:
            line = line.strip().replace('\r', '').replace(' ', '')
            parts = line.split('\t')
            if len(parts) >= 2:
                col_A, col_B = parts[0], parts[1]
                content_headers[col_B] = col_A
except FileNotFoundError:
    print(f"Error: Could not open file {header_list}")
    sys.exit(1)

# Creating dictionary with fasta file content
print("### Obtaining FASTA data")
fasta = {}
try:
    with open(extensive_fasta, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id
            sequence = str(record.seq)
            # Remove trailing asterisk if present
            if sequence.endswith('*'):
                sequence = sequence[:-1]
            fasta[seq_id] = sequence
except FileNotFoundError:
    print(f"Error: Could not open file {extensive_fasta}")
    sys.exit(1)

# Hash for writing new concatenated Fasta
new_fasta_conc = {}

# Main processing loop
for curr_OG in content_OG:
    curr_OG = curr_OG.strip().replace(' ', '')

    print(f"Aligning Block: {curr_OG}...")

    temp_filename = f"{curr_OG}.txt"

    try:
        with open(temp_filename, 'w') as f:
            # Searching inside the selected OGs
            if curr_OG in content_data:
                d = content_data[curr_OG]
                data = d.split(' ')  # All seqs in the OG

                # Focusing on each seq
                for e in data:
                    if e in content_headers:
                        ff = content_headers[e].replace(' ', '')
                        if e in fasta:
                            sqq = fasta[e]
                            f.write(f">{ff}\n{sqq}\n")
    except IOError:
        print(f"Error: Could not create file {temp_filename}")
        continue

    # Creating ALIGNMENT TEMPORAL (MAFFT must be executable in environment)
    aln_filename = f"{curr_OG}.aln"
    try:
        cmd = f"{mafft_exe} --quiet {options} {temp_filename}"
        with open(aln_filename, 'w') as aln_file:
            subprocess.run(cmd, shell=True, stdout=aln_file, check=True)

        # Storing alignment and include it into the current dictionary
        with open(aln_filename, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_id = record.id
                sequence = str(record.seq)

                if seq_id in new_fasta_conc:  # if already present
                    new_fasta_conc[seq_id] += sequence
                else:  # if it is the first
                    new_fasta_conc[seq_id] = sequence

        # Delete temporal files
        os.unlink(temp_filename)
        os.unlink(aln_filename)

    except subprocess.CalledProcessError:
        print(f"Error: MAFFT alignment failed for {curr_OG}")
        # Clean up temp file if it exists
        if os.path.exists(temp_filename):
            os.unlink(temp_filename)
        continue
    except Exception as e:
        print(f"Error processing {curr_OG}: {e}")
        continue

# Final parsing
try:
    with open(OUTPUT, 'w') as out_file:
        print("Assembling...")
        for key, value in new_fasta_conc.items():
            out_file.write(f">{key}\n{value}\n")

    print(f"Done! output available in {OUTPUT}")
except IOError:
    print(f"Error: Could not create output file {OUTPUT}")
    sys.exit(1)

# Function equivalent (though not used in the script)
def uniq(seq):
    seen = set()
    result = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result
