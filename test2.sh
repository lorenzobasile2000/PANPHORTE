#!/bin/bash

# Directory contenente i file GFA
gfa_dir="gfas"

# Directory di output
output_dir="mod"

# Numero massimo di processi in parallelo
num_threads=4

# Creazione della directory di output se non esiste
mkdir -p "$output_dir"

# Funzione per processare un file
process_file() {
    local gfa_file="$1"
    echo "Processing $gfa_file..."
    python3 main.py -i "$gfa_file" -o "$output_dir"
}

# Esportazione della funzione per uso con xargs
export -f process_file
export output_dir

# Trova tutti i file .gfa e processali in parallelo
find "$gfa_dir" -type f -name "*.gfa" | xargs -n 1 -P "$num_threads" bash -c 'process_file "$0"'

echo "Processing complete."
