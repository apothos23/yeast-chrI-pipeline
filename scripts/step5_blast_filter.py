import pandas as pd

blast_file = "results/blast_all_vs_all.tsv"
out_file = "results/step5_homologs.txt"

cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)

print(f"Αρχικές γραμμές BLAST: {len(df)}")

df_diff = df[df["qseqid"] != df["sseqid"]].copy()

# Διαγραφή διπλότυπων ζευγών (A-B και B-A).
# Φτιάχνω μια νέα στήλη 'pair'. Σε αυτή τη στήλη, παίρνω τα δύο ονόματα 
# και τα βάζω πάντα σε ΑΛΦΑΒΗΤΙΚΗ ΣΕΙΡΑ, ενωμένα με μια παύλα (π.χ. YAL001-YAL005).
# Έτσι, είτε διαβάσω το ζευγάρι A-B είτε το B-A, στη στήλη 'pair' θα γραφτεί ακριβώς το ίδιο κείμενο
df_diff["pair"] = df_diff.apply(lambda row: "-".join(sorted([row["qseqid"], row["sseqid"]])), axis=1)

df_unique = df_diff.drop_duplicates(subset=["pair"]).copy()

df_unique = df_unique.drop(columns=["pair"])

df_unique.to_csv(out_file,sep="\t",index=False)

print(f"Τελικά καθαρά και μοναδικά ζεύγη: {len(df_unique)}")
print("\nΔείγμα από τα ζευγάρια που βρέθηκαν:")
print(df_unique[["qseqid", "sseqid", "pident"]].head())