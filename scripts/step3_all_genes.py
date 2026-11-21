import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fasta_file = "data/chrI.fa"
gtf_file = "data/chrI.gtf"
output_file = "results/cds.fa"


col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

df = pd.read_csv(gtf_file, sep="\t", header=None, names=col_names, comment="#")

record = next(SeqIO.parse(fasta_file, "fasta"))
chromosome_name = record.id

cds_all_df = df[df["feature"] == "CDS"].copy()
cds_df = cds_all_df[cds_all_df["seqname"] == chromosome_name].copy()

print(f"Παρέλειψα {len(cds_all_df) - len(cds_df)} γραμμές CDS με διαφορετικό όνομα χρωμοσώματος από το FASTA ({chromosome_name}).")

cds_df["gene_id"] = cds_df["attribute"].str.extract(r'gene_id "([^"]+)"')


cds_df["start"] = cds_df["start"].astype(int)
cds_df["end"] = cds_df["end"].astype(int)

cds_df["length"] = cds_df["end"] - cds_df["start"] + 1

print(cds_df[["seqname", "gene_id", "strand", "start", "end", "length"]].head(10))

chromosome_seq = record.seq.upper()

# Εδώ μαζεύω τα τελικά DNA όλων των γονιδίων
cds_records = []

for gene_id, group in cds_df.groupby("gene_id"):

    # Παίρνω την αλυσίδα (+ ή -) από την πρώτη γραμμή της ομάδας
    strand = group["strand"].iloc[0]
    gene_seq = None

    for index, row in group.iterrows():

        piece = chromosome_seq[row["start"] - 1 : row["end"]]

        if gene_seq is None:
            gene_seq = piece
        else:
            gene_seq += piece

    if strand == "-":
        gene_seq = gene_seq.reverse_complement()

    new_record = SeqRecord(gene_seq, id=gene_id, description="")
    cds_records.append(new_record)

SeqIO.write(cds_records, output_file, "fasta")
print(f"Εξήχθησαν {len(cds_records)} γονίδια στο αρχείο {output_file}.")

final_table = cds_df[["seqname", "gene_id", "strand", "start", "end", "length"]]
output_table_file = "results/step3_genes_table.txt"
final_table.to_csv(output_table_file, sep="\t", index=False)
print(f"Επιτυχία! Ο πίνακας αποθηκεύτηκε στο {output_table_file}.")