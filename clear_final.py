#!/usr/bin/env python3
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pulisce un file GFA eliminando riferimenti a nodi mancanti e correggendo la formattazione delle linee P."
    )
    parser.add_argument("input", help="File GFA in input")
    parser.add_argument("output", help="File GFA in output (pulito)")
    return parser.parse_args()

def extract_node_id(segment):
    """
    Data una stringa segment che rappresenta un nodo,
    se termina con '+' o '-' restituisce solo l'ID numerico.
    In caso contrario, restituisce il segmento così com'è.
    """
    if segment and segment[-1] in ['+', '-']:
        return segment[:-1]
    return segment

def clean_gfa(input_file, output_file):
    valid_nodes = set()
    lines = []

    # Prima lettura: raccogliere tutti gli ID validi dalle linee S
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('S'):
                # Formato atteso: S <node_id> <sequence> [...]
                fields = line.split()
                if len(fields) >= 2:
                    valid_nodes.add(fields[1])
            lines.append(line)

    # Seconda lettura: processare le linee P e L
    output_lines = []
    for line in lines:
        if line.startswith('P'):
            # Formato atteso: P <path_name> <segment_list> <overlap> [opzionali...]
            fields = line.split()
            if len(fields) < 4:
                output_lines.append(line)
                continue

            path_name = fields[1]
            segments = fields[2].split(',')
            overlap = fields[3]

            valid_segments = []
            for seg in segments:
                # Se il segno di orientamento non è presente, aggiungiamo '+'
                if not seg.endswith('+') and not seg.endswith('-'):
                    seg = seg + '+'
                node_id = extract_node_id(seg)
                if node_id in valid_nodes:
                    valid_segments.append(seg)
                # Altrimenti il segmento viene scartato

            if valid_segments:
                new_seg_list = ','.join(valid_segments)
                # Utilizzo del tab come separatore per tutti i campi
                new_fields = [fields[0], path_name, new_seg_list, overlap]
                # Aggiungo eventuali campi opzionali
                if len(fields) > 4:
                    new_fields.extend(fields[4:])
                new_line = "\t".join(new_fields)
                output_lines.append(new_line)
            else:
                # Nessun segmento valido rimasto, quindi si scarta la linea P
                pass

        elif line.startswith('L'):
            # Formato atteso: L <from_node> <from_orient> <to_node> <to_orient> <overlap> [opzionali...]
            fields = line.split()
            if len(fields) < 5:
                output_lines.append(line)
                continue
            from_node = extract_node_id(fields[1])
            to_node = extract_node_id(fields[3])
            if from_node not in valid_nodes or to_node not in valid_nodes:
                # La linea L viene eliminata se almeno uno dei nodi non esiste
                continue
            else:
                output_lines.append(line)
        else:
            # Le altre linee vengono mantenute invariate
            output_lines.append(line)

    # Scriviamo il file di output
    with open(output_file, 'w') as outfile:
        for out_line in output_lines:
            outfile.write(out_line + "\n")

def main():
    args = parse_args()
    clean_gfa(args.input, args.output)

if __name__ == "__main__":
    main()
