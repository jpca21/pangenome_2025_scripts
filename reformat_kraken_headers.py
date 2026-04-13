# !/usr/bin/env python3
"""
reformat_kraken_headers.py
--------------------------
Reformat FASTA sequence headers to Kraken2 custom-database format.

Reads an NCBI assembly_summary.txt to retrieve the taxid associated with
the provided .fna genome file, then emits FASTA with headers formatted as:

    >seqid|kraken:taxid|TAXID [original description]

The Kraken2 standard requires the "kraken:taxid|TAXID" token to either
begin the sequence ID or be immediately preceded by a pipe character (|).

NCBI assembly_summary.txt column layout (tab-delimited, 0-indexed):
    [0]  assembly_accession  (e.g. GCF_965262695.1)
    [5]  taxid               (organism-level NCBI taxonomy ID)
    [6]  species_taxid       (species-level NCBI taxonomy ID)

Usage
-----
    # Single genome  ->  assembly_summary.txt auto-detected in CWD / parent dir
    python reformat_kraken_headers.py GCF_965262695.1_jsCotTube1.Nitrosopumilus_sp_1.1_genomic.fna > output.fna

    # Explicit path to assembly_summary.txt
    python reformat_kraken_headers.py genome.fna /path/to/assembly_summary.txt > output.fna

    # Use species-level taxid instead of organism-level taxid
    python reformat_kraken_headers.py genome.fna --species-taxid > output.fna

    # Batch: process every .fna in the current directory
    for f in *.fna; do
        python reformat_kraken_headers.py "$f" assembly_summary.txt > "kraken_${f}"
    done

Notes
-----
- Biopython is required: pip install biopython
- Accession is extracted from the filename with the pattern GC[AF]_DIGITS.VERSION
- If the versioned accession is absent from the summary, the script retries
  with the unversioned form (e.g. GCF_965262695 without ".1")
- Use --species-taxid to substitute the species-level taxid instead of the
  organism-level taxid (useful when strain-level taxids are too specific)
"""

from __future__ import annotations

import sys
import os
import re
import argparse
from typing import Dict, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_assembly_summary(path):
    # type: (str) -> Dict[str, Dict[str, str]]
    """
    Parse NCBI assembly_summary.txt.

    Returns
    -------
    dict  {assembly_accession: {"taxid": str, "species_taxid": str}}

    The file is tab-delimited; lines beginning with '#' are skipped.
    Relevant 0-indexed columns:
        [0]  assembly_accession
        [5]  taxid
        [6]  species_taxid
    """
    mapping = {}  # type: Dict[str, Dict[str, str]]
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            accession     = cols[0].strip()
            taxid         = cols[5].strip()
            species_taxid = cols[6].strip()
            if accession:
                mapping[accession] = {
                    "taxid":         taxid,
                    "species_taxid": species_taxid,
                }
    return mapping


def get_accession_from_filename(fna_path):
    # type: (str) -> Optional[str]
    """
    Extract GCF/GCA accession (with version) from a .fna filename.

    Example
    -------
    GCF_965262695.1_jsCotTube1.Nitrosopumilus_sp_1.1_genomic.fna
        -> "GCF_965262695.1"
    """
    basename = os.path.basename(fna_path)
    m = re.match(r"(GC[AF]_\d+\.\d+)", basename)
    return m.group(1) if m else None


def find_assembly_summary_auto(fna_path):
    # type: (str) -> Optional[str]
    """
    Search for assembly_summary.txt in the following locations (in order):
        1. Directory that contains the .fna file
        2. One directory level above that
        3. Current working directory
    """
    fna_dir = os.path.dirname(os.path.abspath(fna_path))
    candidates = [
        os.path.join(fna_dir, "assembly_summary.txt"),
        os.path.join(os.path.dirname(fna_dir), "assembly_summary.txt"),
        os.path.join(os.getcwd(), "assembly_summary.txt"),
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return None


def reformat_headers(fna_path, taxid):
    # type: (str, str) -> None
    """
    Stream-parse the FASTA file and write Kraken2-formatted records to stdout.

    Header transformation:
        Before:  >NC_012345.1 Organism name, complete genome
        After:   >NC_012345.1|kraken:taxid|562 Organism name, complete genome

    The original sequence ID is preserved; "|kraken:taxid|TAXID" is appended
    to it so Kraken2 can parse the taxonomy node during database construction.
    Any description text following the ID is kept unchanged.
    """
    for record in SeqIO.parse(fna_path, "fasta"):
        original_id = record.id
        # Everything after the first whitespace-separated token on the header line
        description = record.description[len(original_id):].strip()

        new_id = "{0}|kraken:taxid|{1}".format(original_id, taxid)
        new_record = SeqRecord(
            seq         = record.seq,
            id          = new_id,
            description = description,
        )
        SeqIO.write(new_record, sys.stdout, "fasta")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    # type: () -> argparse.ArgumentParser
    p = argparse.ArgumentParser(
        prog="reformat_kraken_headers.py",
        description=(
            "Reformat FASTA headers for Kraken2 custom database construction "
            "using taxids from an NCBI assembly_summary.txt."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "fna",
        metavar="genome.fna",
        help=(
            "Input genomic FASTA file.  The GCF/GCA accession must be "
            "present in the filename (standard NCBI FTP naming convention)."
        ),
    )
    p.add_argument(
        "summary",
        metavar="assembly_summary.txt",
        nargs="?",
        default=None,
        help=(
            "Path to the NCBI assembly_summary.txt file.  "
            "If omitted the script searches the .fna directory, its parent, "
            "and the current working directory."
        ),
    )
    p.add_argument(
        "--species-taxid",
        action="store_true",
        default=False,
        help=(
            "Use the species-level taxid (column 7, 1-indexed) instead of the "
            "organism-level taxid (column 6, 1-indexed).  Useful when your "
            "taxonomy dump does not include strain-level nodes."
        ),
    )
    return p


def main():
    # type: () -> None
    parser = build_parser()
    args   = parser.parse_args()

    # ---- Validate .fna file ------------------------------------------------
    if not os.path.isfile(args.fna):
        parser.error("FASTA file not found: {0}".format(args.fna))

    # ---- Locate assembly_summary.txt ---------------------------------------
    summary_path = args.summary
    if summary_path is None:
        summary_path = find_assembly_summary_auto(args.fna)
    if summary_path is None or not os.path.isfile(summary_path):
        parser.error(
            "assembly_summary.txt not found.\n"
            "  * Place it in the same directory as the .fna file, or\n"
            "  * Pass its path as the second positional argument."
        )

    # ---- Extract accession from filename -----------------------------------
    accession = get_accession_from_filename(args.fna)
    if not accession:
        parser.error(
            "Could not extract a GCF/GCA accession from filename:\n"
            "  {0}\n"
            "Expected NCBI FTP naming: GCF_XXXXXXXXX.V_<assembly_name>_genomic.fna".format(
                args.fna
            )
        )

    # ---- Load assembly summary ---------------------------------------------
    print("[INFO] Assembly summary  : {0}".format(summary_path), file=sys.stderr)
    acc_map = load_assembly_summary(summary_path)

    # ---- Look up taxid (try versioned, then unversioned) -------------------
    taxid_key = "species_taxid" if args.species_taxid else "taxid"
    entry = acc_map.get(accession)

    if entry is None:
        # Strip version suffix: GCF_965262695.1  -> GCF_965262695
        accession_base = re.sub(r"\.\d+$", "", accession)
        entry = acc_map.get(accession_base)
        if entry:
            print(
                "[WARN] Versioned accession '{0}' not found; "
                "matched unversioned '{1}'.".format(accession, accession_base),
                file=sys.stderr,
            )

    if entry is None:
        sys.exit(
            "[ERROR] Accession '{0}' not found in:\n"
            "  {1}\n"
            "Check that the assembly_summary.txt corresponds to this dataset.".format(
                accession, summary_path
            )
        )

    taxid = entry[taxid_key]
    if not taxid or not taxid.isdigit():
        sys.exit(
            "[ERROR] Invalid taxid value '{0}' for accession '{1}'.\n"
            "        Check column 6 (taxid) in the assembly summary.".format(
                taxid, accession
            )
        )

    print("[INFO] Accession         : {0}".format(accession), file=sys.stderr)
    print("[INFO] TaxID ({0:14s}): {1}".format(taxid_key, taxid), file=sys.stderr)

    # ---- Reformat and stream to stdout -------------------------------------
    reformat_headers(args.fna, taxid)


if __name__ == "__main__":
    main()
