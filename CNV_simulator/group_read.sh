#!/usr/bin/env bash
# Script per elaborare tutti i file FASTA in una cartella di input definita
# e salvare i risultati in una cartella di output definita.
# Le cartelle INPUT_DIR e OUTPUT_DIR sono hardcoded nello script.

set -euo pipefail

# Directory di input (contiene file .fasta o .fa)
INPUT_DIR="t2t/10k/reads_10000"
# Directory di output (saranno creati i file processati)
OUTPUT_DIR="t2t/10k/reads_100000"

# Verifica che la cartella di input esista
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Errore: la cartella di input '$INPUT_DIR' non esiste."
  exit 1
fi

# Crea la cartella di output se non esiste
mkdir -p "$OUTPUT_DIR"

# Loop su tutti i file FASTA (estensioni .fasta e .fa)
for infile in "$INPUT_DIR"/*.fa; do
  if [[ -f "$infile" ]]; then
    base=$(basename "$infile")
    name="${base%.*}"
    outfile="$OUTPUT_DIR/${name}.fa"

    echo "Processing $infile → $outfile"
    python3 group_reads.py -i "$infile" -o "$outfile" -n 10
  fi
done

echo "Tutti i file sono stati elaborati."
