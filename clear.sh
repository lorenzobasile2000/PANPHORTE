#!/bin/bash

# Definizione delle cartelle di input e output
INPUT_DIR="10kmod"
OUTPUT_DIR="10kmodc"
# Numero massimo di processi in parallelo (modifica questo valore in base alla tua CPU)
num_threads=4

# Creazione della cartella di output se non esiste
mkdir -p "$OUTPUT_DIR"

# Funzione per processare un file GFA con lo script Python
process_file() {
    local file="$1"
    local filename
    filename=$(basename -- "$file")
    # Rimuovo l'estensione .gfa
    local filename_noext="${filename%.gfa}"
    # Costruisco il nome del file di output con _c aggiunto prima di .gfa
    local outname="${OUTPUT_DIR}/${filename_noext}_c.gfa"

    echo "Elaborazione di: $filename"
    python3 clear_final.py "$file" "$outname"
    echo "Elaborato: $filename -> $(basename "$outname")"
}

# Esportazione della funzione affinché xargs possa utilizzarla
export -f process_file
export OUTPUT_DIR

# Trova tutti i file .gfa nella cartella INPUT_DIR e processali in parallelo
find "$INPUT_DIR" -type f -name "*.gfa" | xargs -n 1 -P "$num_threads" bash -c 'process_file "$0"'

echo "Processo completato!"
