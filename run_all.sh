#!/bin/bash

set -e

echo "=== Starting the dna2protein pipeline ==="

# 1. Delete the old directory to demonstrate that everything is built from scratch
echo "[1/6] Cleaning and creating the results directory..."
rm -rf results
mkdir results

# 2. Run the Python scripts for statistics and genes
echo "[2/6] Steps 1 & 3: Statistics and Gene Coordinates..."
python scripts/step1_stats.py
python scripts/step3_all_genes.py

# 3. Translation
echo "[3/6] Step 4: Translating genes into proteins..."
python scripts/step4_proteins.py

# 4. BLAST commands (run directly in Bash, as they were run manually)
echo "[4/6] Step 5a: Running BLAST (Database Creation & Search)..."
makeblastdb -in results/proteins.fa -dbtype prot -out results/blastdb
blastp -query results/proteins.fa -db results/blastdb -evalue 1e-5 -outfmt "6 qseqid sseqid pident length evalue bitscore" -out results/blast_all_vs_all.tsv

# 5. BLAST filtering
echo "[5/6] Step 5b: Filtering BLAST results..."
python scripts/step5_blast_filter.py

# 6. Plots
echo "[6/6] Step 6: Creating GC and Length plots..."
python scripts/step6_plots.py

echo "=== The pipeline completed successfully! Check the results/ directory ==="