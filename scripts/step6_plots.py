import matplotlib.pyplot as plt
from Bio import SeqIO

fasta_file = "data/chrI.fa"
proteins_file = "results/proteins.fa"

# --- PLOT 1: GC content per 5,000 bases ---

record = next(SeqIO.parse(fasta_file,"fasta"))
chromosome_seq = record.seq.upper()

window_size = 5000
gc_values = []
positions = []

for i in range(0, len(chromosome_seq), window_size):
    # Exclude the final incomplete window to avoid generating a false peak at the end
    if i + window_size > len(chromosome_seq):
        continue

    window = chromosome_seq[i : i + window_size]
    c_count = window.count("C")
    g_count = window.count("G")

    gc_percent = (c_count + g_count) / len(window) * 100
    gc_values.append(gc_percent)
    positions.append(i)

plt.figure(figsize=(10, 5))
plt.plot(positions, gc_values, color="blue", linewidth=1.5)
plt.title("GC Percentage Across Chromosome I")
plt.xlabel("Genomic Position (bp)")
plt.ylabel("GC Percentage (%)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("results/step6_gc_plot.png")
plt.close()

print("The GC plot was saved to results/step6_gc_plot.png")

# --- PLOT 2: Histogram of protein lengths ---

lengths = []

for record in SeqIO.parse(proteins_file, "fasta"):
    lengths.append(len(record.seq))

plt.figure(figsize=(10, 5))
# The bins=20 option divides the data into 20 columns; edgecolor adds borders to the columns
plt.hist(lengths, bins=20, color="green", edgecolor="black", alpha=0.7)
plt.title("Protein Length Distribution")
plt.xlabel("Protein Length (Amino Acids)")
plt.ylabel("Number of Proteins (Frequency)")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("results/step6_length_histogram.png")
plt.close()

print("The protein-length histogram was saved to results/step6_length_histogram.png")