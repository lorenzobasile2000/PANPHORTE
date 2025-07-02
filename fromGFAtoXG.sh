#!/bin/bash

# Cartelle di input e output
INPUT_DIR="Graph-Morphology-Changer-by-RegexPy/1000modc"
OUTPUT_DIR="Graph-Morphology-Changer-by-RegexPy/1000modcxg"

# Numero massimo di processi paralleli
num_threads=4

# Creazione della cartella di output se non esiste
mkdir -p "$OUTPUT_DIR"

# Funzione per processare un file
process_file() {
    local file="$1"
    local filename
    filename=$(basename -- "$file")
    # Sostituisco l'estensione .gfa con .xg
    local filename_noext="${filename%.gfa}"
    local outname="${OUTPUT_DIR}/${filename_noext}.xg"

    echo "Start conversione: $filename"
    vg convert -g "$file" -x > "$outname"
    echo "Completato: $filename -> $(basename "$outname")"
}

# Esportazione della funzione per l'uso con xargs
export -f process_file
export OUTPUT_DIR

# Trova tutti i file .gfa nella cartella INPUT_DIR e processali in parallelo
find "$INPUT_DIR" -type f -name "*.gfa" | xargs -n 1 -P "$num_threads" bash -c 'process_file "$0"'

echo "Processo di conversione completato!"
