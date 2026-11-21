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
print(f"Αρχικές γραμμές του πίνακα: {len(df)}")

# Φιλτράρισμα: κρατάμε μόνο τις εγγραφές που αντιστοιχούν στο ίδιο χρωμόσωμα με το FASTA
filtered_df = df[df["seqname"] == chromosome_name].copy()
print(f"Παρέλειψα {len(df) - len(filtered_df)} γραμμές με διαφορετικό όνομα χρωμοσώματος από το FASTA ({chromosome_name}).")

df = filtered_df.sort_values(by=["gene_id", "start"])

protein_records = []
stats_list = [] # Εδώ θα μαζέψω τα δεδομένα για τον νέο πίνακα

for gene_id, group in df.groupby("gene_id"):

    full_DNA = Seq("")

    for index, row in group.iterrows():

        piece = chromosome_seq[row["start"]-1 : row["end"]]
        full_DNA += piece

    if group["strand"].iloc[0] == "-":
        
        full_DNA = full_DNA.reverse_complement()

    protein_seq = full_DNA.translate()

    # Αφαιρούμε το τελικό '*' που είναι το stop codon για να έχουμε το σωστό μήκος πρωτεΐνης
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

# 1. Γράφω όλες τις πρωτεΐνες στο αρχείο results/proteins.fa
SeqIO.write(protein_records, out_fasta, "fasta")

# 2. Φτιάχνω τον πίνακα στατιστικών και τον αποθηκεύω
stats_df = pd.DataFrame(stats_list)
stats_df.to_csv(out_table, sep="\t", index=False)

print(f"Επιτυχία! Μεταφράστηκαν {len(protein_records)} γονίδια και αποθηκεύτηκαν στο {out_fasta}")

#Έλεγχος εξαιρέσεων

bad_starts = stats_df[stats_df["starts_with_M"]==False]
bad_stops = stats_df[stats_df["stop_codons"]>0]

print(f"\n---- ΑΝΑΦΟΡΑ ΕΞΑΙΡΕΣΕΩΝ ---")

if not bad_starts.empty:
    print(f"Βρέθηκαν {len(bad_starts)} γονίδια που δεν ξεκινούν με M:")
    print(bad_starts["gene_id"].tolist())
else:
    print("Όλα τα γονίδια ξεκινούν σωστά με Μεθειονίνη (M)!")

if not bad_stops.empty:
    print(f"Βρέθηκαν {len(bad_stops)} γονίδια που περιέχουν Stop Codons (*):")
    print(bad_stops["gene_id"].tolist())
else:
    print("Κανένα γονίδιο δεν έχει πρόωρα stop codons!")