#!/bin/bash

# Definizione delle cartelle di input e output
INPUT_DIR="transfer_dwarf3/CNV_simulator/t2t/1000s/msa"
VG_DIR="transfer_dwarf3/CNV_simulator/t2t/1000s/vgs"
GFA_DIR="transfer_dwarf3/CNV_simulator/t2t/1000s/gfas"

# Creazione delle cartelle di output se non esistono
mkdir -p "$VG_DIR"
mkdir -p "$GFA_DIR"

# Loop su tutti i file .fa nella cartella di input
for file in "$INPUT_DIR"/*.fa; do
    # Estrazione del nome base del file senza estensione
    filename=$(basename -- "$file")
    filename_noext="${filename%.fa}"
    
    # Primo comando: creazione del file .vg
    vg construct -M "$file" -m 60000000 > "$VG_DIR/$filename_noext.vg"
    
    # Secondo comando: conversione del file .vg in .gfa
    vg view "$VG_DIR/$filename_noext.vg" > "$GFA_DIR/$filename_noext.gfa"
    
    echo "Elaborato: $filename"
done

echo "Processo completato!"