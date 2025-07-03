import argparse
import random
import csv

# Parser per gli argomenti della riga di comando
parser = argparse.ArgumentParser(description="Process CNV insertions into a sequence.")
parser.add_argument("-i", "--input", required=True, help="Input file with the sequence")
parser.add_argument("-o", "--output", required=True, help="Output file to generate")
parser.add_argument("-c", "--csv", required=True, help="CSV file with CNV parameters")
args = parser.parse_args()

file_input = args.input
file_output = args.output
csv_file = args.csv

# Funzione per leggere una sequenza dal file
def leggi_sequenza(file_path):
    try:
        with open(file_path, 'r') as file:
            contenuto = file.read().replace('\n', '').strip()  # Rimuove '\n' e spazi vuoti
            # Verifica che il contenuto sia composto solo da caratteri ACGT
            if not all(carattere in "ACGT" for carattere in contenuto):
                raise ValueError("Il file contiene caratteri non validi. Solo 'A', 'C', 'G', 'T' sono ammessi.")
            return contenuto
    except FileNotFoundError:
        print(f"Errore: Il file '{file_path}' non è stato trovato.")
    except ValueError as e:
        print(f"Errore: {e}")

# Funzione per generare una stringa base e stringhe derivate
def genera_stringa_base(lunghezza):
    return ''.join(random.choice("ACGT") for _ in range(lunghezza))

def genera_stringhe(x, y, n_stringhe=8):
    stringa_base = genera_stringa_base(x)
    stringhe_derivate = []
    max_ripetizioni = y // x

    for _ in range(n_stringhe):
        ripetizioni = random.randint(1, max_ripetizioni - 1)
        contenuto = stringa_base * ripetizioni
        padding = '-' * (y - len(contenuto))
        stringhe_derivate.append(contenuto + padding)
    
    stringhe_derivate.append('-' * y)
    max_stringa = stringa_base * max_ripetizioni
    stringhe_derivate.append(max_stringa + '-' * (y - len(max_stringa)))
    
    return stringa_base, stringhe_derivate

# Funzione per scrivere un file con N operazioni di inserimento CNV
def scrivi_file_con_N_CNV(output_path, sequenza, operazioni_CNV):
    try:
        with open(output_path, 'w') as file:
            for i in range(len(operazioni_CNV[0]["stringhe"])):
                seq_modificata = sequenza
                offset_totale = 0
                for op in operazioni_CNV:
                    R = op["R"]
                    stringhe = op["stringhe"]
                    posizione_corrente = R + offset_totale
                    if posizione_corrente < len(seq_modificata):
                        seq_modificata = (
                            seq_modificata[:posizione_corrente]
                            + stringhe[i]
                            + seq_modificata[posizione_corrente:]
                        )
                    else:
                        seq_modificata = seq_modificata + stringhe[i]
                    offset_totale += len(stringhe[i])
                file.write(f">Haplotype:{i + 1}\n")
                file.write(f"{seq_modificata}\n")
        print(f"File generato con successo: {output_path}")
    except Exception as e:
        print(f"Errore durante la scrittura del file: {e}")

# Funzione per leggere il CSV con i parametri di inserimento
def leggi_parametri_csv(csv_file_path):
    operazioni_CNV = []
    try:
        with open(csv_file_path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                R = int(row["R"])
                x = int(row["x"])
                y = int(row["y"])
                _, stringhe_derivate = genera_stringhe(x, y, n_stringhe=8)
                operazioni_CNV.append({
                    "R": R,
                    "stringhe": stringhe_derivate
                })
        return operazioni_CNV
    except FileNotFoundError:
        print(f"Errore: Il file CSV '{csv_file_path}' non è stato trovato.")
    except Exception as e:
        print(f"Errore durante la lettura del file CSV: {e}")
        return []

# Esecuzione principale
sequenza = leggi_sequenza(file_input)
if sequenza:
    operazioni_CNV = leggi_parametri_csv(csv_file)
    if operazioni_CNV:
        scrivi_file_con_N_CNV(file_output, sequenza, operazioni_CNV)
