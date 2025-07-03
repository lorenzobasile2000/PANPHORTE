#!/bin/bash
# Script: run_graphaligner.sh
# Descrizione:
#   - Per ogni file .gfa nella directory specificata da GFA_DIR:
#       * Estrae il nome base (ad es. "chr1" da "chr1.gfa")
#       * Utilizza il corrispondente file di reads nella directory REDAS_DIR (es. "chr1.fa")
#       * Esegue il comando GraphAligner:
#             GraphAligner -g $GFA_DIR/chr1.gfa -f $REDAS_DIR/chr1.fa -a $RESULTS_DIR/chr1_res.gaf -b 35 -x vg -t 6
#       * Dal log generato estrae il numero totale di reads allineate
#   - I risultati vengono salvati in un file CSV (OUTPUT_CSV) contenente una riga per ciascun grafico.
#
# Parametri configurabili:
GFA_DIR="t2t/10k/gfas"           # Directory contenente i file .gfa
REDAS_DIR="t2t/10k/reads_102"       # Directory contenente i file .fa
RESULTS_DIR="t2t/10k/results_pre_10000"   # Directory dove salvare gli output degli allineamenti
OUTPUT_CSV="t2t/10k/results10000_2.csv"  # File CSV di output

# Crea la directory RESULTS se non esiste
[ -d "$RESULTS_DIR" ] || mkdir -p "$RESULTS_DIR"

# Svuota o crea il file CSV di output
> "$OUTPUT_CSV"

# Itera su tutti i file .gfa nella directory GFA_DIR
for graph in "$GFA_DIR"/*.gfa; do
    # Estrai il nome base (ad es. "chr1" da "chr1.gfa")
    base_name=$(basename "$graph" .gfa)
    
    # Definisci il file FASTA corrispondente nella directory REDAS_DIR (ad es. "REDAS/chr1.fa")
    fasta_file="${REDAS_DIR}/${base_name}.fa"
    
    # Verifica che il file FASTA esista
    if [[ ! -f "$fasta_file" ]]; then
        echo "File ${fasta_file} non trovato per ${graph}, salto..."
        continue
    fi
    
    echo "Eseguo GraphAligner per ${base_name}..."
    
    # Definisci il file di output degli allineamenti nella directory RESULTS_DIR
    alignment_output="${RESULTS_DIR}/${base_name}_res.gaf"
    
    # Esegui GraphAligner e cattura l'output (stdout e stderr)
    output=$(GraphAligner -g "$graph" -f "$fasta_file" -a "$alignment_output" -b 35 -x vg -t 7 2>&1)

    # Estrai il numero di input reads; si assume una linea tipo:
    # "Input reads: 10 (900bp)"
    input_reads=$(echo "$output" | grep -i "Input reads:" | sed -E 's/.*Input reads:[[:space:]]*([0-9]+).*/\1/')
    
    # Se il numero non viene estratto, segnala con "N/A"
    if [[ -z "$input_reads" ]]; then
         echo "Non sono riuscito ad estrarre il numero di input reads per ${base_name}"
         input_reads="N/A"
    fi
    
    # Estrai il numero di reads allineate; si assume che l'output contenga una linea tipo:
    # "End-to-end alignments: 12345"
    num_reads=$(echo "$output" | grep -i "End-to-end alignments:" | sed -E 's/.*End-to-end alignments:[[:space:]]*([0-9]+).*/\1/')
    
    # Se il numero non viene estratto, segnala con "N/A"
    if [[ -z "$num_reads" ]]; then
         echo "Non sono riuscito ad estrarre il numero di reads allineate per ${base_name}"
         num_reads="N/A"
    fi
    echo "Elaborato: ${base_name}"
    
    
    # Scrivi il risultato (nome base, input reads, end-to-end alignments) nel file CSV
    echo "${base_name},${input_reads},${num_reads}" >> "$OUTPUT_CSV"
done

echo "Elaborazione completata. I risultati sono stati salvati in ${OUTPUT_CSV}."
