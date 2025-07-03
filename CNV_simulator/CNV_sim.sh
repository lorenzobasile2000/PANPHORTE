#!/bin/bash

# Definizione delle cartelle di input
CHR_DIR="t2t/chr"
CSV_DIR="t2t/10k/csv_tmp"
OUTPUT_DIR="t2t/10k/msa"

# Numero massimo di processi in parallelo
num_threads=2

# Creazione della cartella di output se non esiste
mkdir -p "$OUTPUT_DIR"

# Funzione per eseguire il programma Python
process_file() {
    local base_name="$1"
    local chr_file="$CHR_DIR/$base_name.txt"
    local csv_file="$CSV_DIR/$base_name.csv"
    local output_file="$OUTPUT_DIR/${base_name}_msa_final.fa"
    
    if [[ -f "$chr_file" && -f "$csv_file" ]]; then
        echo "Processing $base_name..."
        python3 CNV_sim.py -i "$chr_file" -o "$output_file" -c "$csv_file"
    else
        echo "Skipping $base_name: missing file(s)"
    fi
}

# Esportazione della funzione per uso con xargs
export -f process_file
export CHR_DIR CSV_DIR OUTPUT_DIR

# Trova i file .txt nella cartella chr e processali in parallelo
find "$CHR_DIR" -type f -name "*.txt" | sed 's#.*/##' | sed 's/\.txt$//' | xargs -n 1 -P "$num_threads" bash -c 'process_file "$0"'

echo "Processing complete."
