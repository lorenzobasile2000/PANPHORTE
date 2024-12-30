import re
import subprocess
import json
import random
import copy
import numpy as np

# Parameters
gfa_file = "MHC_msa_final.gfa"
output_gfa = "out.gfa"
out_json = "OUT_JSON"
fasta = "FASTA"
file_path = "OUT_JSON"
repetition_length = 5


def read_gfa(file_gfa):
    """
    Legge un file GFA e restituisce i record in dizionari strutturati.
    """
    header = []
    nodes = {}
    paths = {}
    links = []

    with open(file_gfa, 'r') as gfa:
        for line in gfa:
            line = line.strip()
            if not line:
                continue

            tipo_record = line[0]

            if tipo_record == 'H':
                header.append(line)

            elif tipo_record == 'S':
                seq = line.split('\t')
                if len(seq) >= 3:
                    id_node = seq[1]
                    sequence = seq[2]
                    nodes[id_node] = sequence

            elif tipo_record == 'P':
                path = line.split('\t')
                if len(path) >= 3:
                    id_path = path[1]
                    segments = [seg[:-1] if seg[-1] in "+-" else seg for seg in path[2].split(',')]
                    paths[id_path] = segments

            elif tipo_record == 'L':
                link = line.split('\t')
                if len(link) >= 5:
                    source = link[1]
                    source_orient = link[2]
                    target = link[3]
                    target_orient = link[4]
                    overlap = link[5]
                    links.append((source, source_orient, target, target_orient, overlap))

    return header, nodes, paths, links


def new_random_id(existing_node_ids, prefix="FUSION"):
    """
    Genera un ID unico per un nuovo nodo, assicurandosi che non sia presente nell'insieme degli ID esistenti.

    :param existing_node_ids: Set di ID già esistenti
    :param prefix: Prefisso per il nuovo ID
    :return: Nuovo ID unico
    """
    while True:
        random_id = f"{prefix}_{random.randint(1, 1000000)}"
        if random_id not in existing_node_ids:
            return random_id


def regex(sequence):
    pattern = re.compile(r'(.+?)\1+')
    matches = pattern.finditer(sequence)

    tandem_repetitions = []
    for match in matches:
        repetition = match.group(1)
        if len(repetition) > repetition_length:
            times = len(match.group(0)) // len(repetition)
            position = match.start()
            tandem_repetitions.append((repetition, times, position))

    # if tandem_repetitions:
    #     for repetition, times, position in tandem_repetitions:
    #         print(f"{repetition} repeated {times} times at position {position}\n")
    # else:
    #     print("No tandem repetition find.\n")

    return tandem_repetitions


header, nodes, paths, links = read_gfa(gfa_file)
# Command to execute Bubblegun
command = ["BubbleGun", "-g", gfa_file, "bchains", "--bubble_json", out_json, "--fasta", fasta]

result = subprocess.run(command, capture_output=True, text=True)

if result.returncode == 0:
    print("Command executed with success!")
    print("Output:", result.stdout)
else:
    print("Error with the execution of the command.")
    print("Error:", result.stderr)

# Opening JSON file
with open(file_path, "r") as file:
    data = json.load(file)

for chains in data:
    chain = data[chains]
    print("Chain ID:", chain["chain_id"])
    print("Chain Ends:", chain["ends"])
    print()
    for bubble in chain["bubbles"]:
        print(f"Bubble ID: {bubble['id']}")
        print(f"Type: {bubble['type']}")
        print(f"Ends: {bubble['ends']}")
        print(f"Inside nodes: {bubble['inside']}")
        print()
        start_node = bubble['ends'][0]
        final_node = bubble['ends'][1]

        if len(bubble['inside']) > 1:  # If not only insertion bubble
            bubble_repetitions = {}
            used = {haplotype: 0 for haplotype in paths}
            haplotypes = []
            # Check if some inside nodes are contiguous, in this case they need to merge into a single pathway
            for inside_node in bubble['inside']:
                for path in paths:
                    if used[path] == 0:
                        if inside_node in paths[path]:
                            contiguous_node = [path, inside_node]
                            for other_node in bubble['inside']:
                                if other_node in paths[path] and other_node != inside_node:
                                    contiguous_node.append(other_node)
                            if len(contiguous_node) > 1:
                                haplotypes.append(contiguous_node)
                                used[path] = 1

            # Find the repetitions
            for haplotype in haplotypes:
                sequence = ''
                if isinstance(haplotype, str):
                    sequence = sequence + nodes[haplotype]
                    # print(f"Calling RegexPy with: {sequence} for node {haplotype}")
                    bubble_repetitions[haplotype] = regex(sequence)
                else:
                    path = haplotype[0]
                    path_of_sequence = []
                    for node in paths[path]:  # TODO: si può semplificare invertendo la logica?
                        if node in haplotype:
                            path_of_sequence.append(node)
                            sequence = sequence + nodes[node]
                    # print(f"Calling RegexPy with: {sequence} \n for nodes {path_of_sequence}")
                    bubble_repetitions[haplotype[0], tuple(path_of_sequence)] = regex(sequence)
            selected_repetition = ()
            original_repetition = copy.deepcopy(bubble_repetitions)
            # Analysis of repetitions
            for node_with_repetition in bubble_repetitions:
                for repetition in bubble_repetitions[node_with_repetition]:
                    if len(repetition) > 0:
                        fusible_nodes = list()
                        fusible_nodes.append(node_with_repetition)
                        for node in bubble_repetitions:
                            if node != node_with_repetition:
                                # Check if any other node contains the same repetition
                                for node_repetition in bubble_repetitions[node]:
                                    if len(node_repetition[0]) > 1:
                                        if node_repetition[0] == repetition[0] and node != node_with_repetition:
                                            fusible_nodes.append(node)
                                            bubble_repetitions[node].remove(node_repetition)
                                # Or if any other node contains the repeated sequence
                                if node not in fusible_nodes and isinstance(node, str):
                                    if repetition[0] in nodes[node]:
                                        fusible_nodes.append(node)
                                elif node not in fusible_nodes:
                                    sequence = ''
                                    for inside_node in node[1]:
                                        sequence = sequence + nodes[inside_node]
                                    if repetition[0] in sequence:
                                        fusible_nodes.append(node)
                        if len(fusible_nodes) > 1:
                            if len(selected_repetition) == 0 or len(repetition[0]) > len(selected_repetition[0][0]):
                                selected_repetition = (repetition, fusible_nodes)
                            # else:
                            #     print(
                            #         f"No structure modification proposed for repetition {repetition} of node {node_with_repetition}")
                            #     # bubble_repetitions[node_with_repetition].remove(repetition)
            # Modify the morphology by reconstructing the original sequences of the different pathway and comparing it
            print(f"{selected_repetition[1]} can be fused due to {selected_repetition[0]}!")
            # Reconstruct original sequences
            fusible_sequences = {}
            for fusible_path in selected_repetition[1]:
                sequence = ''
                for node in fusible_path[1]:
                    sequence += nodes[node]
                fusible_sequences[fusible_path[0]] = sequence
            # New node for the repetition, valid for each haplotype
            new_rep_id = new_random_id(nodes, "REP")
            nodes[new_rep_id] = selected_repetition[0][0]
            up_flk_dict = {}
            dw_flk_dict = {}
            for haplotype in selected_repetition[1]:
                current_sequence = fusible_sequences[haplotype[0]]
                # Select the right repetition
                current_rep = tuple()
                for rep in original_repetition[haplotype]:
                    if rep[0] == selected_repetition[0][0]:
                        current_rep = rep
                # It means that the node doesn't contain repetition itself, but the repetition is present one time
                if len(current_rep) < 2:
                    current_rep = (selected_repetition[0][0], 1, current_sequence.find(selected_repetition[0][0]))
                # New node for upstream and downstream flanking sequence, if needed
                # If repetition doesn't start at position 0
                if current_rep[2] > 0:
                    # This is the upstream flanking sequence
                    up_flk_seq = current_sequence[:current_rep[2]]
                    # If I didn't create a node for this upstream sequence yet
                    if up_flk_seq not in up_flk_dict:
                        # In this case I need to create a new node, link them to the starting node of the
                        # bubble and to the repeated node
                        new_up_flk_id = new_random_id(nodes, "UP_FLK")
                        up_flk_dict[up_flk_seq] = new_up_flk_id
                        nodes[new_up_flk_id] = up_flk_seq
                        links.append((start_node, '+', new_up_flk_id, '+', '0M'))
                        links.append((new_up_flk_id, '+', new_rep_id, '+', '0M'))
                    else:
                        new_up_flk_id = up_flk_dict[up_flk_seq]
                else:
                    # In this case I need to link starting node to repeated node
                    new_up_flk_id = 'NONE'
                    links.append((start_node, '+', new_rep_id, '+', '0M'))
                # This is the downstream flanking sequence
                # rep[0] = sequence of the repetition
                # rep[1] = times
                # rep[2] = starting position
                tmp = len(current_sequence)
                dw_flk_position = current_rep[2] + len(current_rep[0]) * current_rep[1]
                # If the repetition doesn't end at the sequence length
                if dw_flk_position < len(current_sequence):
                    dw_flk_seq = current_sequence[dw_flk_position:]
                    # If I didn't create a node for this downstream sequence yet
                    if dw_flk_seq not in dw_flk_dict:
                        # In this case I need to create a new node, link them to the repeated node and to the
                        # final node of the bubble
                        new_dw_flk_id = new_random_id(nodes, "DW_FLK")
                        dw_flk_dict[dw_flk_seq] = new_dw_flk_id
                        nodes[new_dw_flk_id] = dw_flk_seq
                        links.append((new_rep_id, '+', new_dw_flk_id, '+', '0M'))
                        links.append((new_dw_flk_id, '+', final_node, '+', '0M'))
                    else:
                        new_dw_flk_id = dw_flk_dict[dw_flk_seq]
                else:
                    # In this case I need to link repeated node to the final node
                    new_dw_flk_id = 'NONE'
                    links.append((new_rep_id, '+', final_node, '+', '0M'))
                # Now I'm ready to modify the pathway
                current_path = haplotype[0]
                # Find the position of the starting node
                index = paths[current_path].index(start_node)
                # If upstream flanking node
                if new_up_flk_id != 'NONE':
                    paths[current_path].insert(index+1, new_up_flk_id)
                    index += 1
                # Adding repeated node for each time the repetition is present
                times = current_rep[1]
                for i in range(1, times):
                    paths[current_path].insert(index + i, new_rep_id)
                # If downstream flanking node
                if new_dw_flk_id != 'NONE':
                    paths[current_path].insert(index + 1, new_dw_flk_id)
            # Remove old nodes
            deleted_nodes = []
            for path in selected_repetition[1]:
                for node in path[1]:
                    if node in nodes:
                        del nodes[node]
                        deleted_nodes.append(node)
            # Cleaning up old links
            for link in links:
                if link[0] in deleted_nodes or link[2] in deleted_nodes:
                    links.remove(link)
        else:
            print("No further analysis needed!")
        print("-----")

with open(output_gfa, 'w') as gfa:
    # Scrivi l'header
    for header_line in header:
        gfa.write(f"{header_line}\n")

    # Scrivi i nodi
    for node_id, sequence in nodes.items():
        gfa.write(f"S\t{node_id}\t{sequence}\n")

    # Scrivi i percorsi
    for path_id, segments in paths.items():
        gfa.write(f"P\t{path_id}\t{','.join(segments)}\t*\n")

    # Scrivi i collegamenti
    for link in links:
        gfa.write(f"L\t{link[0]}\t{link[1]}\t{link[2]}\t{link[3]}\t{link[4]}\n")

print(f"New GFA written to {output_gfa}")
