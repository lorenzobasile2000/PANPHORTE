import re
from collections import defaultdict, deque


class Link:
    def __init__(self, start, end, id):
        self.start = start
        self.end = end
        self.id = id


def leggi_gfa(file_gfa):
    links = []
    nodes = {}
    count = 0

    with open(file_gfa, 'r') as file:
        for linea in file:
            linea = linea.strip()
            if not linea:
                continue

            tipo_record = linea[0]

            if tipo_record == 'L':
                dati = linea.split('\t')
                if len(dati) == 6:
                    start = dati[1]
                    end = dati[3]
                    link = Link(start, end, count)
                    links.append(link)
                    count += 1
                else:
                    print(f"Riga non valida: {linea}")

            elif tipo_record == 'S':
                dati = linea.split('\t')
                if len(dati) >= 3:
                    id_node = dati[1]
                    sequence = dati[2]
                    nodes[id_node] = sequence  # Aggiungi il nodo e la sua sequenza

    return links, nodes


def trova_bolle(link_list):
    grafo = defaultdict(list)
    for link in link_list:
        grafo[link.start].append(link.end)

    bolle = []

    def trova_percorsi(nodo_inizio):
        percorsi = defaultdict(list)
        coda = deque([[nodo_inizio]])

        while coda:
            percorso = coda.popleft()
            ultimo_nodo = percorso[-1]

            if ultimo_nodo in grafo:
                for next_nodo in grafo[ultimo_nodo]:
                    nuovo_percorso = percorso + [next_nodo]
                    coda.append(nuovo_percorso)
                    percorsi[next_nodo].append(nuovo_percorso)

        return percorsi

    for nodo in grafo:
        if len(grafo[nodo]) > 1:
            percorsi_da_nodo = trova_percorsi(nodo)
            nodi_finali = defaultdict(list)
            for finale, percorsi in percorsi_da_nodo.items():
                if len(percorsi) > 1:  # Modifica qui: richiediamo almeno 3 percorsi
                    nodi_finali[finale].extend(percorsi)

            for nodo_finale, percorsi in nodi_finali.items():
                if len(percorsi) > 1:  # Assicuriamoci di avere almeno 3 percorsi per considerare una bolla
                    bolle.append((nodo, nodo_finale, percorsi))
                    break    #Se troviamo una bolla, smettiamo di cercare altri percorsi per questo nodo

    return bolle


def regex(sequenza_genoma):
    pattern = re.compile(r'(.+?)\1+')
    matches = pattern.finditer(sequenza_genoma)

    tandem_repetitions = []
    for match in matches:
        ripetizione = match.group(1)
        numero_di_volte = len(match.group(0)) // len(ripetizione)
        posizione = match.start()
        tandem_repetitions.append((ripetizione, numero_di_volte, posizione))

    if tandem_repetitions:
        for ripetizione, numero_di_volte, posizione in tandem_repetitions:
            print(f"{ripetizione} repeated {numero_di_volte} times at position {posizione}\n")
    else:
        print("No tandem repetition find.\n")

    return tandem_repetitions


# Esegui la funzione con un esempio di file GFA
file_gfa = "example3.gfa"
links, nodes = leggi_gfa(file_gfa)

# Ordina la lista di links in base al campo "start"
links_ordinati = sorted(links, key=lambda link: link.start)

# Identificare le bolle
bolle_trovate = trova_bolle(links_ordinati)

# Stampare le bolle trovate con le sequenze concatenate
if bolle_trovate:
    for bolla in bolle_trovate:
        nodo_inizio, nodo_fine, percorsi = bolla
        print(f"Bolla trovata tra {nodo_inizio} e {nodo_fine}:")
        for percorso in percorsi:
            if len(percorso) > 2:  # Assicurarsi che ci siano nodi intermedi
                # Escludi il primo (nodo_inizio) e l'ultimo (nodo_fine) nodo
                sequenza = ''.join(nodes[nodo] for nodo in percorso[1:-1])
                if sequenza:  # Controllo aggiuntivo per evitare stringhe vuote
                    print(f"Sequenza esclusi nodi estremi: {sequenza}")
                    regex(sequenza)
            else:
                print("Percorso troppo corto per escludere nodi estremi.")
else:
    print("Nessuna bolla trovata.")
