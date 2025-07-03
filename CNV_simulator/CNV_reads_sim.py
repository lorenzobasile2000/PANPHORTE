import argparse
import random
import csv

def read_fasta(file_path):
    """
    Legge un file FASTA e restituisce un dizionario in cui 
    le chiavi sono gli header (senza '>') e i valori le sequenze.
    """
    haplotypes = {}
    with open(file_path, 'r') as f:
        header = None
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    haplotypes[header] = seq
                header = line[1:]
                seq = ""
            else:
                seq += line
        if header is not None:
            haplotypes[header] = seq
    return haplotypes

def read_csv(cnv_csv):
    """
    Legge il file CSV e restituisce una lista di operazioni CNV.
    Ogni operazione è un dizionario contenente:
      - R: posizione nella sequenza originale
      - x: lunghezza dell'elemento ripetitivo (unità)
      - y: lunghezza totale dell'inserto (in caso di massima ripetizione)
    """
    operations = []
    with open(cnv_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                R = int(row["R"])
                x = int(row["x"])
                y = int(row["y"])
                operations.append({"R": R, "x": x, "y": y})
            except Exception as e:
                print("Errore nella lettura del CSV:", e)
    return operations

def select_max_haplotype(haplotypes):
    """
    Seleziona il haplotype con la sequenza più lunga,
    che assumiamo essere quello generato con il massimo inserto.
    """
    max_header = max(haplotypes, key=lambda k: len(haplotypes[k]))
    return max_header, haplotypes[max_header]

def recover_repetitive_elements(max_seq, operations):
    """
    Dato il haplotype massimo e la lista di operazioni, estrae per ogni operazione
    l'elemento ripetitivo: si preleva dalla porzione inserita (di lunghezza y)
    i primi x caratteri, che corrispondono alla unità di ripetizione.
    """
    rep_elements = []
    offset = 0
    for op in operations:
        pos = op["R"] + offset
        inserted_segment = max_seq[pos: pos + op["y"]]
        rep_unit = inserted_segment[:op["x"]]
        rep_elements.append(rep_unit)
        offset += op["y"]
    return rep_elements

def recover_original_sequence(max_seq, operations):
    """
    Ricostruisce la sequenza originale a partire dal haplotype massimo,
    rimuovendo gli inserimenti. Si assume che le operazioni siano ordinate
    in base alla posizione R crescente.
    """
    original = ""
    offset = 0
    last_idx = 0
    for op in operations:
        pos = op["R"] + offset
        # Aggiunge la porzione della sequenza fino all'inserto
        original += max_seq[last_idx:pos]
        # Aggiorna l'offset e il punto di inizio per la prossima iterazione
        offset += op["y"]
        last_idx = pos + op["y"]
    original += max_seq[last_idx:]
    return original

def generate_new_haplotypes(original_seq, operations, rep_elements, n_haplotypes=20):
    """
    A partire dalla sequenza originale, reinserisce per ciascuna operazione
    l'elemento ripetitivo recuperato, ripetendolo un numero casuale di volte.
    Il numero massimo di ripetizioni è dato da (y // x).
    """
    haplotypes = []
    for i in range(n_haplotypes):
        seq_mod = original_seq
        offset = 0
        for op, rep_unit in zip(operations, rep_elements):
            max_copies = op["y"] // op["x"]
            copies = random.randint(1, max_copies)
            inserted = rep_unit * copies
            pos = op["R"] + offset
            if pos < len(seq_mod):
                seq_mod = seq_mod[:pos] + inserted + seq_mod[pos:]
            else:
                seq_mod = seq_mod + inserted
            offset += len(inserted)
        haplotypes.append(seq_mod)
    return haplotypes

def write_fasta_reads(output_path, haplotypes, max_read_length=10000):
    """
    Scrive le sequenze degli haplotype nel file FASTA di output, 
    suddividendole in reads di massimo max_read_length caratteri.
    Ogni read avrà header nel formato: >R:<numero_haplotype>:<numero_read>
    """
    with open(output_path, "w") as f:
        for i, hap_seq in enumerate(haplotypes):
            # Divide la sequenza in pezzi di lunghezza massima max_read_length
            for j in range(0, len(hap_seq), max_read_length):
                read = hap_seq[j:j+max_read_length]
                read_number = j // max_read_length + 1
                f.write(f">R:{i+1}:{read_number}\n")
                f.write(read + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Recupera elemento ripetitivo da un file FASTA (10 haplotype) e genera 20 nuovi haplotype "
                    "inserendo casualmente un numero di ripetizioni dell'elemento (senza padding). "
                    "Le sequenze finali vengono suddivise in reads di max 1000 caratteri."
    )
    parser.add_argument("-f", "--fasta", required=True, help="File FASTA di input (con 10 haplotype)")
    parser.add_argument("-c", "--csv", required=True, help="File CSV con parametri CNV (R, x, y)")
    parser.add_argument("-o", "--output", required=True, help="File di output FASTA con 20 nuovi haplotype")
    args = parser.parse_args()
    
    # Lettura del file FASTA
    haplotypes = read_fasta(args.fasta)
    if not haplotypes:
        print("Nessun haplotype trovato nel file FASTA.")
        return
    
    # Selezione del haplotype con massima lunghezza (assunto essere quello con inserimenti massimi)
    max_header, max_seq = select_max_haplotype(haplotypes)
    print(f"Utilizzo il haplotype '{max_header}' per recuperare gli elementi ripetitivi.")
    
    # Lettura del file CSV
    operations = read_csv(args.csv)
    if not operations:
        print("Nessuna operazione trovata nel CSV. Terminazione.")
        return
    
    # Recupero degli elementi ripetitivi (unità) dal haplotype massimo
    rep_elements = recover_repetitive_elements(max_seq, operations)
    
    # Ricostruzione della sequenza originale rimuovendo gli inserimenti
    original_seq = recover_original_sequence(max_seq, operations)
    
    # Generazione di 20 nuovi haplotype a partire dalla sequenza originale,
    # inserendo per ciascuna operazione il ripetitivo (unità) ripetuto casualmente
    new_haplotypes = generate_new_haplotypes(original_seq, operations, rep_elements, n_haplotypes=20)
    
    # Scrittura dell'output in un file FASTA, con suddivisione in reads da max 1000 char
    write_fasta_reads(args.output, new_haplotypes)
    print(f"File FASTA con 20 nuovi haplotype generato con successo: {args.output}")

if __name__ == "__main__":
    main()
