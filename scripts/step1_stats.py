from Bio import SeqIO

fasta_file = "data/chrI.fa"

output_file = "results/step1_stats.txt"

for record in SeqIO.parse(fasta_file, "fasta-blast"):
    dna = record.seq
    count_A = dna.count("A")
    count_C = dna.count("C")
    count_G = dna.count("G")
    count_T = dna.count("T")
    total_length = len(dna)
    gc_percent = (count_G+count_C)/(count_A+count_C+count_G+count_T)*100
    print(f"Chromosome length: {total_length} bases")
    print(f"A: {count_A}, C: {count_C}, G: {count_G}, T: {count_T}")
    print(f"GC percentage: {gc_percent}%")

report = (
    f"Total length: {total_length} bases/n"
    f"A: {count_A}, C: {count_C}, T: {count_T}, G: {count_G}/n"
    f"GC percentage: {gc_percent:.2f}%/n"
)

with open(output_file, "w", encoding="utf-8") as f:
    f.write(report)