##############################################################################
##############################################################################
#
#  USAGE: python add_peptide_position_ultra_fast_v2.py -f DB.fasta -r diann.report.tsv
#
#  ULTRA-OPTIMIZED VERSION with ROBUST FASTA ID parsing and parallel processing
#
#  LOW-MEMORY MODE: the report is streamed and processed in row batches
#  (see --batch-size) and the output is written incrementally, so peak
#  memory no longer scales with the full report size. On a very
#  memory-constrained machine, use a small --batch-size together with
#  --jobs 1 (no worker pool at all -> only a single copy of the FASTA
#  lookup tables in memory).
#
##############################################################################
##############################################################################


import os
import sys
from Bio import SeqIO
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
from tqdm import tqdm
import argparse
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser(
        description = "Adding protein position to reported peptides - ULTRA OPTIMIZED")



###  INPUT handling

parser.add_argument('-f', "--fasta", help = "Path to fasta-file used as DIANN library", required = True)
parser.add_argument('-r', "--report", help = "DIANN output report,  report.tsv or report.parquet allowed", required = True)
parser.add_argument('-j', "--jobs", help = "Number of parallel jobs (default: use all cores). Use 1 to disable multiprocessing entirely and keep only a single copy of the FASTA lookup tables in memory.", type=int, default=None)
parser.add_argument('-b', "--batch-size", dest="batch_size", help = "Number of report rows streamed/processed at a time (default: 200000). Lower this on low-memory machines when processing large (~1GB+) report files.", type=int, default=200000)

args = parser.parse_args()

if args.batch_size < 1:
    parser.error("--batch-size must be >= 1")

# Determine number of cores
n_jobs = args.jobs if args.jobs else max(1, cpu_count() - 1)
print(f"Using {n_jobs} parallel job(s), batch size {args.batch_size} rows", file=sys.stderr)

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


# Populated once per worker process via the Pool initializer (or once in the
# main process for -j 1) instead of being re-pickled and re-sent for every
# chunk, which previously multiplied both memory and IPC overhead.
_protein_id_to_seq = None
_fasta_id_variants = None


def init_worker(protein_dict, variants_dict):
    global _protein_id_to_seq, _fasta_id_variants
    _protein_id_to_seq = protein_dict
    _fasta_id_variants = variants_dict


def process_row_chunk(row_dicts):
    """
    Process a chunk of rows (in a worker process, or serially in the main
    process when -j 1 is used). Relies on _protein_id_to_seq /
    _fasta_id_variants having been set via init_worker().
    """
    protein_id_to_seq = _protein_id_to_seq
    fasta_id_variants = _fasta_id_variants

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


def process_batch(row_dicts, pool):
    """
    Split one in-memory batch of rows into chunks and process them, either
    across the worker pool or serially in the main process (pool is None).
    """
    if pool is not None:
        chunk_size = max(1, len(row_dicts) // (n_jobs * 4))
        chunks = [row_dicts[i:i + chunk_size] for i in range(0, len(row_dicts), chunk_size)]
        results = pool.map(process_row_chunk, chunks)
    else:
        results = [process_row_chunk(row_dicts)]

    expanded_data = []
    not_found_ids = set()
    for chunk_result, chunk_not_found in results:
        expanded_data.extend(chunk_result)
        not_found_ids.update(chunk_not_found)

    return expanded_data, not_found_ids


# process diann report in either parquet or tsv format, streaming it in
# row batches and writing output incrementally so memory use stays bounded
# by --batch-size rather than by the full report size

diann_report = args.report
diann_file_name = os.path.split(args.report)[1]

pool = Pool(n_jobs, initializer=init_worker, initargs=(protein_id_to_seq, fasta_id_variants)) if n_jobs > 1 else None
if pool is None:
    init_worker(protein_id_to_seq, fasta_id_variants)

total_rows = 0
total_expanded = 0
all_not_found = set()

try:
    if diann_report.endswith(".parquet"):

        print("\nDetected .parquet file\n")

        parquet_file = pq.ParquetFile(diann_report)
        source_total_rows = parquet_file.metadata.num_rows
        n_batches_est = max(1, -(-source_total_rows // args.batch_size))
        print(f"Report contains {source_total_rows} rows, streaming in batches of {args.batch_size}", file=sys.stderr)

        output_file = diann_file_name.replace(".parquet", ".pos.parquet")
        writer = None
        output_schema = None

        for batch in tqdm(parquet_file.iter_batches(batch_size=args.batch_size),
                           total=n_batches_est, desc="batches processed", ascii=True):
            batch_df = batch.to_pandas()
            total_rows += batch_df.shape[0]
            row_dicts = batch_df.to_dict('records')
            del batch_df

            expanded_data, not_found_ids = process_batch(row_dicts, pool)
            all_not_found.update(not_found_ids)
            del row_dicts

            out_df = pd.DataFrame(expanded_data)
            del expanded_data
            out_df["prot_ID"] = out_df["prot_ID"].astype("string")
            out_df["pept_start"] = out_df["pept_start"].astype("Int64")
            out_df["pept_stop"] = out_df["pept_stop"].astype("Int64")

            table = pa.Table.from_pandas(out_df, preserve_index=False)
            del out_df

            if writer is None:
                output_schema = table.schema
                writer = pq.ParquetWriter(output_file, output_schema)
            else:
                table = table.cast(output_schema)

            writer.write_table(table)
            total_expanded += table.num_rows
            del table

        if writer is not None:
            writer.close()

        print(f"Original rows: {total_rows}")
        print(f"Expanded rows: {total_expanded}")
        if all_not_found:
            print(f"Warning: {len(all_not_found)} unique protein IDs not found in FASTA", file=sys.stderr)
            print(f"First 10 missing IDs: {list(all_not_found)[:10]}", file=sys.stderr)
        print("Output written to:", output_file)

    elif diann_report.endswith(".tsv"):

        print("\nDetected .tsv file\n")
        print(f"Streaming TSV file in batches of {args.batch_size} rows...", file=sys.stderr)

        new_tsv_file_name = diann_file_name.replace(".tsv", ".pos.tsv")
        first_batch = True

        reader = pd.read_csv(diann_report, sep='\t', low_memory=False, chunksize=args.batch_size)
        for batch_df in tqdm(reader, desc="batches processed", ascii=True):
            total_rows += batch_df.shape[0]
            row_dicts = batch_df.to_dict('records')
            del batch_df

            expanded_data, not_found_ids = process_batch(row_dicts, pool)
            all_not_found.update(not_found_ids)
            del row_dicts

            out_df = pd.DataFrame(expanded_data)
            del expanded_data
            out_df.to_csv(new_tsv_file_name, sep='\t', index=False,
                           mode='w' if first_batch else 'a', header=first_batch)
            total_expanded += out_df.shape[0]
            del out_df
            first_batch = False

        print(f"Original rows: {total_rows}")
        print(f"Expanded rows: {total_expanded}")
        if all_not_found:
            print(f"Warning: {len(all_not_found)} unique protein IDs not found in FASTA", file=sys.stderr)
            print(f"First 10 missing IDs: {list(all_not_found)[:10]}", file=sys.stderr)
        print("Output written to:", new_tsv_file_name)

    else:
        parser.error("--report must end with .tsv or .parquet")

finally:
    if pool is not None:
        pool.close()
        pool.join()
