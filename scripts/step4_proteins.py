import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta_file = "data/chrI.fa"
table_file = "results/step3_genes_table.txt"
out_fasta = "results/proteins.fa"
out_table = "results/step4_protein_stats.txt"

record = next(SeqIO.parse(fasta_file,"fasta"))
chromosome_name = record.id
chromosome_seq = record.seq.upper()

df = pd.read_csv(table_file,sep="\t")
print(f"Initial table rows: {len(df)}")

# Filter to retain only records that correspond to the same chromosome as the FASTA
filtered_df = df[df["seqname"] == chromosome_name].copy()
print(f"Skipped {len(df) - len(filtered_df)} rows with a chromosome name different from the FASTA ({chromosome_name}).")

df = filtered_df.sort_values(by=["gene_id", "start"])

protein_records = []
stats_list = [] # Collect data for the new table here

for gene_id, group in df.groupby("gene_id"):

    full_DNA = Seq("")

    for index, row in group.iterrows():

        piece = chromosome_seq[row["start"]-1 : row["end"]]
        full_DNA += piece

    if group["strand"].iloc[0] == "-":
        
        full_DNA = full_DNA.reverse_complement()

    protein_seq = full_DNA.translate()

    # Remove the terminal '*' stop codon to obtain the correct protein length
    prot_str = str(protein_seq)
    if prot_str.endswith("*"):
        prot_str = prot_str[:-1]

    protein_seq = Seq(prot_str)

    new_record = SeqRecord(protein_seq, id=gene_id, description="")
    
    protein_records.append(new_record)

    length = len(prot_str)
    starts_with_M = prot_str.startswith("M")
    stop_codons = prot_str.count("*")

    stats_list.append({
        "gene_id": gene_id,
        "length": length,
        "starts_with_M": starts_with_M,
        "stop_codons": stop_codons
    })

# 1. Write all proteins to results/proteins.fa
SeqIO.write(protein_records, out_fasta, "fasta")

# 2. Create the statistics table and save it
stats_df = pd.DataFrame(stats_list)
stats_df.to_csv(out_table, sep="\t", index=False)

print(f"Success! Translated {len(protein_records)} genes and saved them to {out_fasta}")

#Exception check

bad_starts = stats_df[stats_df["starts_with_M"]==False]
bad_stops = stats_df[stats_df["stop_codons"]>0]

print(f"\n---- EXCEPTION REPORT ---")

if not bad_starts.empty:
    print(f"Found {len(bad_starts)} genes that do not start with M:")
    print(bad_starts["gene_id"].tolist())
else:
    print("All genes start correctly with Methionine (M)!")

if not bad_stops.empty:
    print(f"Found {len(bad_stops)} genes that contain Stop Codons (*):")
    print(bad_stops["gene_id"].tolist())
else:
    print("No gene contains premature stop codons!")