#!/usr/bin/env python3
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Pulisce un file GFA eliminando riferimenti a nodi mancanti e correggendo la separazione dei campi nelle linee P.")
    parser.add_argument("input", help="File GFA in input")
    parser.add_argument("output", help="File GFA in output (pulito)")
    return parser.parse_args()

def extract_node_id(segment):
    """
    Dato un segmento che può essere in forma 'ID+' o 'ID-', ritorna solo l'ID.
    Se non termina con + o -, restituisce il segmento così com'è.
    """
    if segment and segment[-1] in ['+', '-']:
        return segment[:-1]
    return segment

def clean_gfa(input_file, output_file):
    valid_nodes = set()
    lines = []

    # Prima lettura: raccogliere tutti gli ID dei nodi dalle linee S
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('S'):
                # Consideriamo il formato: S <node_id> <sequence> [...]
                fields = line.rstrip().split()
                if len(fields) >= 2:
                    valid_nodes.add(fields[1])
            lines.append(line.rstrip())

    # Seconda lettura: processare le linee P e L
    output_lines = []
    for line in lines:
        if line.startswith('P'):
            # Struttura attesa: P <path_name> <segment_list> <overlap> [opzionali...]
            fields = line.split()
            if len(fields) < 4:
                output_lines.append(line)
                continue

            path_name = fields[1]
            segments = fields[2].split(',')
            overlap = fields[3]

            # Manteniamo solo i segmenti validi
            valid_segments = []
            for seg in segments:
                node_id = extract_node_id(seg)
                if node_id in valid_nodes:
                    valid_segments.append(seg)
                else:
                    # Segmento rimosso (nodo non presente in S)
                    pass

            # Se esistono segmenti validi, compongo la nuova linea P con campi separati da tab
            if valid_segments:
                new_seg_list = ','.join(valid_segments)
                new_fields = [fields[0], path_name, new_seg_list, overlap]
                # Aggiungo eventuali campi opzionali, mantenendo la separazione a tab
                if len(fields) > 4:
                    new_fields.extend(fields[4:])
                new_line = "\t".join(new_fields)
                output_lines.append(new_line)
            else:
                # Se non rimane nessun segmento valido, si scarta la linea P
                pass

        elif line.startswith('L'):
            # Struttura attesa: L <from_node> <from_orient> <to_node> <to_orient> <overlap> [opzionali...]
            fields = line.split()
            if len(fields) < 5:
                output_lines.append(line)
                continue
            from_node = extract_node_id(fields[1])
            to_node = extract_node_id(fields[3])
            if from_node not in valid_nodes or to_node not in valid_nodes:
                # Elimino la linea se uno dei due nodi non è presente
                continue
            else:
                output_lines.append(line)
        else:
            # Le altre linee (S, header, ecc.) vengono lasciate invariate
            output_lines.append(line)

    # Scrittura del file di output
    with open(output_file, 'w') as outfile:
        for out_line in output_lines:
            outfile.write(out_line + "\n")

def main():
    args = parse_args()
    clean_gfa(args.input, args.output)

if __name__ == "__main__":
    main()
