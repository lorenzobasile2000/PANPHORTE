def process_gfa(input_path, output_path):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('P\tHaplotype:'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    nodes = parts[2].split(',')
                    # Aggiungi '+' se non è già presente
                    updated_nodes = [node if node.endswith('+') else f"{node}+" for node in nodes]
                    parts[2] = ','.join(updated_nodes)
                    new_line = '\t'.join(parts)
                    outfile.write(new_line + '\n')
                else:
                    outfile.write(line)  # se non rispetta il formato atteso
            else:
                outfile.write(line)

# Esempio di utilizzo
process_gfa('100mod/chrY_mod.gfa', 'output.gfa')
