#!/bin/bash

echo "=== Ξεκινάει το Pipeline dna2protein ==="

# 1. Διαγραφή του παλιού φακέλου (για να αποδείξουμε ότι όλα χτίζονται από την αρχή)
echo "[1/6] Καθαρισμός και δημιουργία φακέλου results..."
rm -rf results
mkdir results

# 2. Εκτέλεση των Python scripts για στατιστικά και γονίδια
echo "[2/6] Βήμα 1 & 3: Στατιστικά και Συντεταγμένες Γονιδίων..."
python scripts/step1_stats.py
python scripts/step3_all_genes.py

# 3. Μετάφραση
echo "[3/6] Βήμα 4: Μετάφραση γονιδίων σε πρωτεΐνες..."
python scripts/step4_proteins.py

# 4. Εντολές BLAST (κατευθείαν στο bash, όπως τις έτρεξες χειροκίνητα)
echo "[4/6] Βήμα 5α: Εκτέλεση BLAST (Δημιουργία Βάσης & Αναζήτηση)..."
makeblastdb -in results/proteins.fa -dbtype prot -out results/blastdb
blastp -query results/proteins.fa -db results/blastdb -evalue 1e-5 -outfmt "6 qseqid sseqid pident length evalue bitscore" -out results/blast_all_vs_all.tsv

# 5. Φιλτράρισμα BLAST
echo "[5/6] Βήμα 5β: Φιλτράρισμα αποτελεσμάτων BLAST..."
python scripts/step5_blast_filter.py

# 6. Γραφήματα
echo "[6/6] Βήμα 6: Δημιουργία Γραφημάτων GC και Μηκών..."
python scripts/step6_plots.py

echo "=== Το Pipeline ολοκληρώθηκε επιτυχώς! Έλεγξε τον φάκελο results/ ==="