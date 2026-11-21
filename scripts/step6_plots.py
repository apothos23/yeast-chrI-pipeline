import matplotlib.pyplot as plt
from Bio import SeqIO

fasta_file = "data/chrI.fa"
proteins_file = "results/proteins.fa"

# --- ΓΡΑΦΗΜΑ 1: GC Content ανά 5.000 βάσεις ---

record = next(SeqIO.parse(fasta_file,"fasta"))
chromosome_seq = record.seq.upper()

window_size = 5000
gc_values = []
positions = []

for i in range(0, len(chromosome_seq), window_size):
    # Αφήνουμε έξω το τελευταίο, ατελές παράθυρο για να μην δημιουργείται ψεύτικο peak στο τέλος
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
plt.title("Ποσοστό GC κατά μήκος του Χρωμοσώματος I")
plt.xlabel("Γενομική Θέση (bp)")
plt.ylabel("Ποσοστό GC (%)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("results/step6_gc_plot.png")
plt.close()

print("Το γράφημα GC αποθηκεύτηκε στο results/step6_gc_plot.png")

# --- ΓΡΑΦΗΜΑ 2: Ιστόγραμμα μηκών πρωτεϊνών ---

lengths = []

for record in SeqIO.parse(proteins_file, "fasta"):
    lengths.append(len(record.seq))

plt.figure(figsize=(10, 5))
# Το bins=20 σπάει τα δεδομένα σε 20 κολόνες, το edgecolor βάζει περίγραμμα στις κολόνες
plt.hist(lengths, bins=20, color="green", edgecolor="black", alpha=0.7)
plt.title("Κατανομή Μηκών Πρωτεϊνών")
plt.xlabel("Μήκος Πρωτεΐνης (Αμινοξέα)")
plt.ylabel("Αριθμός Γονιδίων (Συχνότητα)")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("results/step6_length_histogram.png")
plt.close()

print("Το ιστόγραμμα μηκών αποθηκεύτηκε στο results/step6_length_histogram.png")