#!/home/momi/mambaforge/bin/python3.10
##############################################################################
##############################################################################
#
#  USAGE: python add_peptide_position_ultra_fast_v2.py DB.fasta diann.report.tsv
#
#  ULTRA-OPTIMIZED VERSION with ROBUST FASTA ID parsing and parallel processing
#
##############################################################################
##############################################################################


import os
import codecs
import sys
from Bio import SeqIO
import pyarrow.parquet as pq
import pandas as pd
from tqdm import tqdm
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial

parser = argparse.ArgumentParser(
        description = "Adding protein position to reported peptides - ULTRA OPTIMIZED")



###  INPUT handling

parser.add_argument('-f', "--fasta", help = "Path to fasta-file used as DIANN library", required = True)
parser.add_argument('-r', "--report", help = "DIANN output report,  report.tsv or report.parquet allowed", required = True)
parser.add_argument('-j', "--jobs", help = "Number of parallel jobs (default: use all cores)", type=int, default=None)


args = parser.parse_args()

# Determine number of cores
n_jobs = args.jobs if args.jobs else max(1, cpu_count() - 1)
print(f"Using {n_jobs} parallel jobs", file=sys.stderr)

fasta_lib = args.fasta


def parse_fasta_id(description):
    """
    Robustly parse FASTA IDs from various formats:
    - >sp|P12345|PROTEIN_NAME
    - >prf|Ricin_E|Ricin_E
    - >tr|A0A0H3GKP4|A0A0H3GKP4_HUMAN
    - >P12345
    - >PROTEIN_NAME

    Returns the protein ID (tries to get the most specific identifier)
    """
    # Remove '>' if present
    desc = description.lstrip('>')

    # Split by pipe
    parts = desc.split('|')

    if len(parts) >= 2:
        # For formats like sp|P12345|NAME or prf|Ricin_E|Ricin_E
        # Take the second part (index 1) which is usually the primary ID
        return parts[1].strip()
    elif len(parts) == 1:
        # For simple format like >P12345 or >PROTEIN_NAME
        # Take the first part and clean it (remove any spaces)
        return parts[0].split()[0].strip()
    else:
        # Fallback: return the whole description
        return desc.strip()


# load fasta to memory - CREATE INDEXED DICTIONARY FOR O(1) LOOKUP
print("Loading FASTA file...", file=sys.stderr)
protein_id_to_seq = {}  # Direct protein_id -> sequence mapping
fasta_id_variants = {}  # Maps all possible ID variants to the canonical ID

fasta_sequences = SeqIO.parse(open(fasta_lib),'fasta')
for fasta in fasta_sequences:
        id, name, descri = fasta.id, fasta.name, fasta.description
        sequence = str(fasta.seq)

        # Use robust parsing
        diann_prot_id = parse_fasta_id(descri)

        protein_id_to_seq[diann_prot_id] = sequence

        # Also store the full ID and name as variants (for matching flexibility)
        fasta_id_variants[id] = diann_prot_id
        fasta_id_variants[name] = diann_prot_id

        # Store all parts of pipe-separated IDs
        if '|' in descri:
            for part in descri.split('|'):
                part = part.strip().lstrip('>')
                if part:
                    fasta_id_variants[part] = diann_prot_id

print(f"Loaded {len(protein_id_to_seq)} protein sequences", file=sys.stderr)
print(f"Created {len(fasta_id_variants)} ID variants for matching", file=sys.stderr)


def process_row_chunk(row_dicts, protein_id_to_seq, fasta_id_variants):
    """
    Process a chunk of rows in parallel
    """
    expanded_data = []
    not_found_ids = set()

    for row_dict in row_dicts:
        stripped_seq = row_dict.get("Stripped.Sequence", "")
        protein_ids_str = row_dict.get("Protein.Ids", "")

        if not protein_ids_str or not stripped_seq:
            # Keep the original row even when it can't be processed
            new_row = row_dict.copy()
            new_row["prot_ID"] = None
            new_row["pept_start"] = None
            new_row["pept_stop"] = None
            expanded_data.append(new_row)
            continue

        protein_ids = protein_ids_str.split(";")
        row_matched = False

        # Process each protein ID
        for prot_id in protein_ids:
            prot_id = prot_id.strip()

            # Strategy 1: Direct lookup
            canonical_id = None
            prot_seq = None

            if prot_id in protein_id_to_seq:
                canonical_id = prot_id
                prot_seq = protein_id_to_seq[prot_id]
            # Strategy 2: Try variant lookup
            elif prot_id in fasta_id_variants:
                canonical_id = fasta_id_variants[prot_id]
                if canonical_id in protein_id_to_seq:
                    prot_seq = protein_id_to_seq[canonical_id]

            if prot_seq is not None:
                pept_start = prot_seq.find(stripped_seq)

                if pept_start != -1:
                    # Create new row
                    new_row = row_dict.copy()
                    new_row["prot_ID"] = canonical_id
                    new_row["pept_start"] = pept_start
                    new_row["pept_stop"] = pept_start + len(stripped_seq)
                    expanded_data.append(new_row)
                    row_matched = True
                else:
                    not_found_ids.add(prot_id)
            else:
                not_found_ids.add(prot_id)

        if not row_matched:
            # No protein ID matched/positioned - keep the original row anyway
            new_row = row_dict.copy()
            new_row["prot_ID"] = None
            new_row["pept_start"] = None
            new_row["pept_stop"] = None
            expanded_data.append(new_row)

    return expanded_data, not_found_ids


# process diann report in either parquet or tsv format

diann_report = args.report
diann_file_name = os.path.split(args.report)[1]


if diann_report.endswith(".parquet"):

    print("\nDetected .parquet file\n")
    dia_df = pd.read_parquet(diann_report, engine='pyarrow')

    print(f"Processing {dia_df.shape[0]} rows in report")

    # Convert to list of dictionaries
    row_dicts = dia_df.to_dict('records')

    # Split into chunks for parallel processing
    chunk_size = max(1, len(row_dicts) // (n_jobs * 4))  # Create more chunks than cores for better load balancing
    chunks = [row_dicts[i:i + chunk_size] for i in range(0, len(row_dicts), chunk_size)]

    print(f"Processing {len(chunks)} chunks in parallel...", file=sys.stderr)

    # Process chunks in parallel
    process_func = partial(process_row_chunk,
                           protein_id_to_seq=protein_id_to_seq,
                           fasta_id_variants=fasta_id_variants)

    with Pool(n_jobs) as pool:
        results = list(tqdm(
            pool.imap(process_func, chunks),
            total=len(chunks),
            desc="chunks processed",
            ascii=True
        ))

    # Flatten results
    expanded_data = []
    all_not_found = set()
    for chunk_result, not_found_ids in results:
        expanded_data.extend(chunk_result)
        all_not_found.update(not_found_ids)

    print(f"Original rows: {dia_df.shape[0]}")
    print(f"Expanded rows: {len(expanded_data)}")
    if all_not_found:
        print(f"Warning: {len(all_not_found)} unique protein IDs not found in FASTA", file=sys.stderr)
        print(f"First 10 missing IDs: {list(all_not_found)[:10]}", file=sys.stderr)

    # Create DataFrame once from all data
    dia_df_expanded = pd.DataFrame(expanded_data)

    #write pandas df back to parquet
    output_file = diann_file_name.replace(".parquet", ".pos.parquet")
    dia_df_expanded.to_parquet(output_file, index=False)
    print("Output written to:", output_file)


elif diann_report.endswith(".tsv"):
    print("\nDetected .tsv file\n")

    # For TSV, use pandas for faster processing
    print("Reading TSV file...", file=sys.stderr)
    dia_df = pd.read_csv(diann_report, sep='\t', low_memory=False)

    print(f"Processing {dia_df.shape[0]} rows in report")

    # Convert to list of dictionaries
    row_dicts = dia_df.to_dict('records')

    # Split into chunks for parallel processing
    chunk_size = max(1, len(row_dicts) // (n_jobs * 4))
    chunks = [row_dicts[i:i + chunk_size] for i in range(0, len(row_dicts), chunk_size)]

    print(f"Processing {len(chunks)} chunks in parallel...", file=sys.stderr)

    # Process chunks in parallel
    process_func = partial(process_row_chunk,
                           protein_id_to_seq=protein_id_to_seq,
                           fasta_id_variants=fasta_id_variants)

    with Pool(n_jobs) as pool:
        results = list(tqdm(
            pool.imap(process_func, chunks),
            total=len(chunks),
            desc="chunks processed",
            ascii=True
        ))

    # Flatten results
    expanded_data = []
    all_not_found = set()
    for chunk_result, not_found_ids in results:
        expanded_data.extend(chunk_result)
        all_not_found.update(not_found_ids)

    print(f"Original rows: {dia_df.shape[0]}")
    print(f"Expanded rows: {len(expanded_data)}")
    if all_not_found:
        print(f"Warning: {len(all_not_found)} unique protein IDs not found in FASTA", file=sys.stderr)
        print(f"First 10 missing IDs: {list(all_not_found)[:10]}", file=sys.stderr)

    # Create DataFrame and write
    dia_df_expanded = pd.DataFrame(expanded_data)

    new_tsv_file_name = diann_file_name.replace(".tsv", ".pos.tsv")
    print("Writing output file...", file=sys.stderr)
    dia_df_expanded.to_csv(new_tsv_file_name, sep='\t', index=False)
    print("Output written to:", new_tsv_file_name)
