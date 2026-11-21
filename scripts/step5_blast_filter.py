import pandas as pd

blast_file = "results/blast_all_vs_all.tsv"
out_file = "results/step5_homologs.txt"

cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)

print(f"Initial BLAST rows: {len(df)}")

df_diff = df[df["qseqid"] != df["sseqid"]].copy()

# Remove duplicate pairs (A-B and B-A).
# Create a new 'pair' column. In this column, take the two names
# and always sort them ALPHABETICALLY, joining them with a dash (e.g., YAL001-YAL005).
# Thus, whether the pair is read as A-B or B-A, the 'pair' column contains exactly the same text
df_diff["pair"] = df_diff.apply(lambda row: "-".join(sorted([row["qseqid"], row["sseqid"]])), axis=1)

df_unique = df_diff.drop_duplicates(subset=["pair"]).copy()

df_unique = df_unique.drop(columns=["pair"])

df_unique.to_csv(out_file,sep="\t",index=False)

print(f"Final clean, unique pairs: {len(df_unique)}")
print("\nSample of pairs found:")
print(df_unique[["qseqid", "sseqid", "pident"]].head())