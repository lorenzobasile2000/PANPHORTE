#!/bin/bash

# Cartelle di input e output
INPUT_DIR="t2t/10k/msa"
CSV_DIR="t2t/10k/csv"
OUTPUT_DIR="t2t/10k/reads_10000"
PYTHON_SCRIPT="CNV_reads_sim.py"  # Modifica con il nome (e path) corretto del tuo script Python

# Numero massimo di processi paralleli
num_threads=1

# Creazione della cartella di output se non esiste
mkdir -p "$OUTPUT_DIR"

# Funzione per processare un file FASTA e il corrispondente CSV
process_file() {
    local fasta_file="$1"
    local filename
    filename=$(basename -- "$fasta_file")
    # Rimuovo l'estensione (es. .fa) per ottenere il nome base
    local filename_noext="${filename%.*}"
    local csv_file="${CSV_DIR}/${filename_noext}.csv"
    local outname="${OUTPUT_DIR}/${filename_noext}.fa"

    # Verifica l'esistenza del file CSV associato
    if [[ ! -f "$csv_file" ]]; then
        echo "Attenzione: file CSV non trovato per $filename. Atteso: ${csv_file}"
        return
    fi

    echo "Avvio elaborazione: $filename"
    python "$PYTHON_SCRIPT" -f "$fasta_file" -c "$csv_file" -o "$outname"
    echo "Completato: $filename -> $(basename "$outname")"
}

# Esportazione della funzione e delle variabili per l'uso con xargs
export -f process_file
export CSV_DIR
export OUTPUT_DIR
export PYTHON_SCRIPT

# Trova tutti i file .fa nella cartella INPUT_DIR ed esegue in parallelo process_file
find "$INPUT_DIR" -type f -name "*.fa" | xargs -n 1 -P "$num_threads" bash -c 'process_file "$0"'

echo "Elaborazione completata!"
