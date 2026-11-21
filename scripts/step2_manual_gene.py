from Bio import SeqIO

fasta_file = "data/chrI.fa"

record = next(SeqIO.parse(fasta_file,"fasta"))
chromosome_seq = record.seq.upper()

#gene1 "YAL067W-A" (+)

start_gtf = 2480
end_gtf = 2704

cds_seq_plus = chromosome_seq[start_gtf-1 : end_gtf]

protein_plus = cds_seq_plus.translate()

print("=== Gene YAL067W-A (+ strand) ===")
print(f"DNA sequence: {cds_seq_plus}")
print(f"Protein: {protein_plus}")

#gene2 "YAL037C-A" (−)

start_gtf_minus = 73429	
end_gtf_minus = 73518

cds_seq_plus_minus = chromosome_seq[start_gtf_minus - 1 : end_gtf_minus]
reversed_cds_seq_plus_minus = cds_seq_plus_minus.reverse_complement()

protein_plus_minus = reversed_cds_seq_plus_minus.translate()

print("=== Gene YAL037C-A (- strand) ===")
print(f"DNA sequence: {reversed_cds_seq_plus_minus}")
print(f"Protein: {protein_plus_minus}")