#!/home/momi/mambaforge/bin/python3.10
##############################################################################
##############################################################################
#
#  USAGE: python add_peptide_position.py DB.fasta diann.report.tsv
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

parser = argparse.ArgumentParser(
        description = "Adding protein position to reported peptides")



###  INPUT handling

parser.add_argument('-f', "--fasta", help = "Path to fasta-file used as DIANN library", required = True)
parser.add_argument('-r', "--report", help = "DIANN output report,  report.tsv or report.parquet allowed", required = True)


args = parser.parse_args()



fasta_lib = args.fasta




# load fasta to memory
fasta_dict = {}

fasta_sequences = SeqIO.parse(open(fasta_lib),'fasta')
for fasta in fasta_sequences:
        id, name, descri,  sequence = fasta.id, fasta.name, fasta.description,  str(fasta.seq)
        diann_prot_id = descri.split("|")[1]
        fasta_dict[descri] = [sequence, diann_prot_id]

for key, value in fasta_dict.items():
    pass
    #print(key, value, file = sys.stderr)




# process diann report in either parquet or tsv format

diann_report = args.report
diann_file_name = os.path.split(args.report)[1]

# match protein.ids with fasta headers
diann_protein_ids = []
out_list = []


if diann_report.endswith(".parquet"):

    pcount = 0
    print("\nDetected .parquet file\n")
    dia_df = pd.read_parquet(diann_report, engine='pyarrow')
    
    print("Processing ", dia_df.shape[0], "found rows in report")
    
    # Create a list to store expanded rows
    expanded_rows = []
    
    for index, row in tqdm(dia_df.iterrows(), ascii = True, desc = "rows processed: ", total=dia_df.shape[0]):
        
        #match peptides to fasta - iterate over each protein ID
        for e in row["Protein.Ids"].split(";"):
            e = e.strip()  # Remove any whitespace
            
            for prot_seq, prot_id in fasta_dict.values():
                if e == prot_id:
                    # Create a copy of the row for each protein ID
                    new_row = row.copy()
                    new_row["prot_ID"] = e
                    strip_peptide_len = len(row["Stripped.Sequence"])
                    pept_start = int(prot_seq.find(row["Stripped.Sequence"]))
                    new_row["pept_start"] = pept_start
                    new_row["pept_stop"] = pept_start + strip_peptide_len
                    
                    expanded_rows.append(new_row)
                    break

    print("Original rows:", dia_df.shape[0])
    
    # Create new dataframe from expanded rows
    dia_df_expanded = pd.DataFrame(expanded_rows)
    
    print("Expanded rows:", dia_df_expanded.shape[0])

    #write pandas df back to parquet pos.parquet
    dia_df_expanded.to_parquet(diann_file_name.replace(".parquet", ".pos.parquet"))
    print("Output written to: ", diann_file_name.replace(".parquet", ".pos.parquet"))


elif diann_report.endswith(".tsv"):
    print("\nDetected .tsv file\n")

    with open(diann_report) as diann_rep:
            first_line = diann_rep.readline().strip().split("\t")
            first_line.append("prot_ID")
            first_line.append("pept_start")
            first_line.append("pept_stop")
            out_list.append(first_line)
            line_count = 1
            #print("\t".join(first_line))
            for line in diann_rep:
                k = 0
                line_count += 1

                line = line.strip().split("\t")
                
                # Split protein IDs and create one output line per protein
                for e in line[3].split(";"):
                    e = e.strip()  # Remove any whitespace
                    print(e, file = sys.stderr)
                    #searching
                    for values in fasta_dict.values():
                        if e in values:
                            #get peptide
                            out_line = line.copy()

                            strip_peptide = line[14]
                            strip_peptide_len = len(line[14])
                            first_hit_start = int(values[0].find(strip_peptide))
                            out_line.append(e)
                            out_line.append(str(first_hit_start))
                            out_line.append(str((first_hit_start + strip_peptide_len)))
                            #print("\t".join(out_line))
                            out_list.append(out_line)
                            k += 1
                            break  # Only take first match per protein ID
                            #print(isinstance(strip_peptide_len, int))
                            #print(strip_peptide)
                            #print(values[0][first_hit_start:(first_hit_start + strip_peptide_len)])
                if k == 0:
                    print("No match found for line:", line, file=sys.stderr)

            print("Total lines processed:", line_count, file=sys.stderr)
            print("Total output lines:", len(out_list), file=sys.stderr)


         # write output
    diann_file_name = os.path.split(diann_report)[1]
    new_tsv_file_name = diann_file_name.replace(".tsv", ".pos.tsv")
    file2 = codecs.open(new_tsv_file_name, "w", "utf-8")
    file2.write('\n'.join("\t".join(i) for i in out_list))
    file2.close()
    print("Output written to:", new_tsv_file_name)

