#!/usr/bin/env python3

def aggiungi_plus_nodi(input_gfa, output_gfa):
    with open(input_gfa, 'r') as fin, open(output_gfa, 'w') as fout:
        for line in fin:
            # Verifica se la riga inizia con 'P '
            if line.startswith('P'):
                # Esempio di riga: "P Haplotype:1 2,3,4,5"
                # Suddividiamo in massimo 3 parti (col 'split(None, 2)'):
                #   1) "P"
                #   2) "Haplotype:1"
                #   3) "2,3,4,5"
                parti = line.strip().split(None, 2)

                # Se non ci sono almeno 3 parti, scriviamo la riga così com'è
                if len(parti) < 3:
                    fout.write(line)
                    continue

                # parti[0] = "P"
                # parti[1] = "Haplotype:<numero>"
                # parti[2] = "lista_di_nodi_separati_da_virgole"
                nodi_str = parti[2]

                # Suddividiamo i nodi separati da virgola
                nodi = nodi_str.split(',')

                # Aggiungiamo un '+' alla fine di ogni nodo, ignorando eventuali stringhe vuote
                nodi_modificati = [n + '+' for n in nodi if n]

                # Ricostruiamo la linea
                nuova_linea = f"{parti[0]} {parti[1]} {','.join(nodi_modificati)}\n"
                fout.write(nuova_linea)
            else:
                # Se la riga non inizia con 'P', la scriviamo senza modifiche
                fout.write(line)


if __name__ == "__main__":
    # Esempio di utilizzo
    input_file = "input.gfa"
    output_file = "output.gfa"
    aggiungi_plus_nodi(input_file, output_file)
    print(f"File trasformato salvato in: {output_file}")
