import argparse
import re
import subprocess
import json
import random
import copy
import os

# -------------------- CLI --------------------
parser = argparse.ArgumentParser(description="Process a GFA file and modify it (W-lines only).")
parser.add_argument("-i", "--input", required=True, help="Input GFA file (GFA v1.1+ with W-lines)")
parser.add_argument("-o", "--output_dir", default=".", help="Output directory for the modified GFA file")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)
gfa_file = args.input
output_gfa = os.path.join(args.output_dir, os.path.splitext(os.path.basename(gfa_file))[0] + "_mod.gfa")

out_json = "OUT_JSON"
fasta = "FASTA"
file_path = "OUT_JSON"
repetition_length = 1

# -------------------- Helpers W-lines --------------------
WALK_TOKEN_RE = re.compile(r'([><])([^><]+)')  # -> ('>', 'segId') or ('<', 'segId')

def parse_walk_string(walk_str):
    """
    Parse W field like '>s11<s12>s13' into a list of (segId, orient) where orient in {'+','-'}.
    """
    result = []
    for m in WALK_TOKEN_RE.finditer(walk_str.strip()):
        arrow, seg = m.groups()
        orient = '+' if arrow == '>' else '-'
        result.append((seg, orient))
    return result

def format_walk_segments(segments):
    """
    Serialize list of (segId, orient) back to W walk string.
    '+' -> '>', '-' -> '<'
    """
    parts = []
    for seg, orient in segments:
        arrow = '>' if orient == '+' else '<'
        parts.append(f"{arrow}{seg}")
    return ''.join(parts)

def find_index_by_seg(segments, seg_id):
    """
    Return first index where segments[idx][0] == seg_id, else raise ValueError.
    """
    for i, (sid, _) in enumerate(segments):
        if sid == seg_id:
            return i
    raise ValueError(f"{seg_id} not found in walk")

def contains_seg(segments, seg_id):
    return any(sid == seg_id for sid, _ in segments)

# -------------------- I/O GFA --------------------
def read_gfa(file_gfa):
    """
    Lettura GFA: header, nodi (S), walk (W), link (L).
    Restituisce:
      header: [str]
      nodes: {nodeId: sequence}
      walks: {walk_id: {'sample':..., 'hap':int, 'seqid':..., 'start':str, 'end':str, 'segments':[(seg,orient), ...]}}
      links: [(from,forient,to,torient,overlap)]
    """
    header = []
    nodes = {}
    walks = {}
    links = []

    walk_auto_id = 0

    with open(file_gfa, 'r') as gfa:
        for line in gfa:
            line = line.rstrip('\n')
            if not line:
                continue
            t = line[0]

            if t == 'H':
                header.append(line)

            elif t == 'S':
                cols = line.split('\t')
                if len(cols) >= 3:
                    node_id = cols[1]
                    sequence = cols[2]
                    nodes[node_id] = sequence

            elif t == 'L':
                cols = line.split('\t')
                if len(cols) >= 6:
                    links.append((cols[1], cols[2], cols[3], cols[4], cols[5]))

            elif t == 'W':
                # W  SampleId  HapIndex  SeqId  SeqStart  SeqEnd  Walk
                cols = line.split('\t')
                if len(cols) >= 7:
                    sample = cols[1]
                    hap = cols[2]
                    seqid = cols[3]
                    seqstart = cols[4]
                    seqend = cols[5]
                    walk_field = cols[6]
                    segments = parse_walk_string(walk_field)
                    walk_id = f"W_{walk_auto_id}:{sample}|{hap}|{seqid}|{seqstart}|{seqend}"
                    walk_auto_id += 1
                    walks[walk_id] = {
                        'sample': sample,
                        'hap': int(hap),
                        'seqid': seqid,
                        'start': seqstart,
                        'end': seqend,
                        'segments': segments
                    }

            # ignore P-lines entirely

    return header, nodes, walks, links

def new_random_id(existing_node_ids, prefix="FUSION"):
    while True:
        rid = f"{prefix}_{random.randint(1, 1_000_000)}"
        if rid not in existing_node_ids:
            return rid

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
    return tandem_repetitions

# -------------------- Main --------------------
header, nodes, walks, links = read_gfa(gfa_file)

# Esegue BubbleGun per le superbubble
command = ["BubbleGun", "-g", gfa_file, "bchains", "--bubble_json", out_json, "--fasta", fasta]
result = subprocess.run(command, capture_output=True, text=True)
if result.returncode == 0:
    print("Command executed with success!")
    print("Output:", result.stdout)
else:
    print("Error with the execution of the command.")
    print("Error:", result.stderr)

# Carica JSON delle bubble
with open(file_path, "r") as file:
    data = json.load(file)

# ---- Corpo logico: come prima, ma agendo su 'walks' ----
innn=0
for chains in data:
    chain = data[chains]
    for bubble in chain["bubbles"]:
        print(str(innn)+"\n")
        innn=innn+1
        start_node = bubble['ends'][0]
        final_node = bubble['ends'][1]
        ok = 0
        # Identifica ordine corretto degli end guardando la prima walk che li contiene
        for w_id in walks:
            for (sid, _) in walks[w_id]['segments']:
                if sid == start_node:
                    ok = 1
                    break
                if sid == final_node:
                    final_node, start_node = start_node, final_node
                    ok = 1
                    break
            if ok:
                break

        if len(bubble['inside']) > 1:
            bubble_repetitions = {}
            used = {w_id: 0 for w_id in walks}
            haplotypes = []

            # Raggruppa nodi interni contigui per walk
            for inside_node in bubble['inside']:
                for w_id, w in walks.items():
                    if used[w_id] == 0 and contains_seg(w['segments'], inside_node):
                        contiguous = [w_id, inside_node]
                        for other_node in bubble['inside']:
                            if other_node != inside_node and contains_seg(w['segments'], other_node):
                                contiguous.append(other_node)
                        if len(contiguous) > 1:
                            haplotypes.append(contiguous)
                            used[w_id] = 1

            # Trova ripetizioni per ciascun haplotype
            for haplotype in haplotypes:
                sequence = ''
                if isinstance(haplotype, str):
                    sequence += nodes[haplotype]
                    bubble_repetitions[haplotype] = regex(sequence)
                else:
                    w_id = haplotype[0]
                    path_nodes = []
                    for (sid, _) in walks[w_id]['segments']:
                        if sid in haplotype:
                            path_nodes.append(sid)
                            sequence += nodes[sid]
                    bubble_repetitions[(w_id, tuple(path_nodes))] = regex(sequence)

            selected_repetition = ()
            original_repetition = copy.deepcopy(bubble_repetitions)

            # Analizza ripetizioni e seleziona il motif migliore da fondere
            for node_with_rep in list(bubble_repetitions.keys()):
                for repetition in list(bubble_repetitions[node_with_rep]):
                    if len(repetition) > 0:
                        fusible_nodes = [node_with_rep]
                        for node in list(bubble_repetitions.keys()):
                            if node != node_with_rep:
                                for node_rep in list(bubble_repetitions[node]):
                                    if len(node_rep[0]) > 1 and node_rep[0] == repetition[0]:
                                        fusible_nodes.append(node)
                                        bubble_repetitions[node].remove(node_rep)
                                if node not in fusible_nodes and isinstance(node, str):
                                    if repetition[0] in nodes[node]:
                                        fusible_nodes.append(node)
                                elif node not in fusible_nodes:
                                    seq2 = ''.join(nodes[n] for n in node[1])
                                    if repetition[0] in seq2:
                                        fusible_nodes.append(node)

                        if len(fusible_nodes) > 1:
                            if (len(selected_repetition) == 0 or
                                len(repetition[0]) > len(selected_repetition[0][0])):
                                selected_repetition = (repetition, fusible_nodes)

            # Modifica struttura della bubble usando il nodo ciclico
            if len(selected_repetition) == 2:
                fusible_sequences = {}
                for fusible in selected_repetition[1]:
                    seq = ''.join(nodes[n] for n in fusible[1])
                    fusible_sequences[fusible[0]] = seq

                # Crea nodo ripetitivo + loop su se stesso
                new_rep_id = new_random_id(nodes, "REP")
                nodes[new_rep_id] = selected_repetition[0][0]
                links.append((new_rep_id, '+', new_rep_id, '+', '0M'))

                up_flk_dict = {}
                dw_flk_dict = {}
                link_between_start_rep = 0
                link_between_rep_end = 0

                for hap in selected_repetition[1]:
                    current_sequence = fusible_sequences[hap[0]]
                    current_rep = tuple()
                    for rep in original_repetition[hap]:
                        if rep[0] == selected_repetition[0][0]:
                            current_rep = rep
                    if len(current_rep) < 2:
                        current_rep = (selected_repetition[0][0], 1, current_sequence.find(selected_repetition[0][0]))

                    # Upstream flanking
                    if current_rep[2] > 0:
                        up_seq = current_sequence[:current_rep[2]]
                        if up_seq not in up_flk_dict:
                            new_up_id = new_random_id(nodes, "UP_FLK")
                            up_flk_dict[up_seq] = new_up_id
                            nodes[new_up_id] = up_seq
                            links.append((start_node, '+', new_up_id, '+', '0M'))
                            links.append((new_up_id, '+', new_rep_id, '+', '0M'))
                        else:
                            new_up_id = up_flk_dict[up_seq]
                    else:
                        new_up_id = 'NONE'
                        if link_between_start_rep == 0:
                            links.append((start_node, '+', new_rep_id, '+', '0M'))
                            link_between_start_rep = 1

                    # Downstream flanking
                    dw_pos = current_rep[2] + len(current_rep[0]) * current_rep[1]
                    if dw_pos < len(current_sequence):
                        dw_seq = current_sequence[dw_pos:]
                        if dw_seq not in dw_flk_dict:
                            new_dw_id = new_random_id(nodes, "DW_FLK")
                            dw_flk_dict[dw_seq] = new_dw_id
                            nodes[new_dw_id] = dw_seq
                            links.append((new_rep_id, '+', new_dw_id, '+', '0M'))
                            links.append((new_dw_id, '+', final_node, '+', '0M'))
                        else:
                            new_dw_id = dw_flk_dict[dw_seq]
                    else:
                        new_dw_id = 'NONE'
                        if link_between_rep_end == 0:
                            links.append((new_rep_id, '+', final_node, '+', '0M'))
                            link_between_rep_end = 1

                    # --- Aggiorna la WALK corrispondente ---
                    current_wid = hap[0]
                    segs = walks[current_wid]['segments']
                    # posizione dello start_node nella walk
                    idx = find_index_by_seg(segs, start_node)

                    # inserisci upstream flanking se creato
                    if new_up_id != 'NONE':
                        segs.insert(idx + 1, (new_up_id, '+'))
                        idx += 1

                    # inserisci il nodo ripetizione 'times' volte
                    times = current_rep[1]
                    insert_at = idx + 1
                    for i in range(times):
                        segs.insert(insert_at + i, (new_rep_id, '+'))
                        idx += 1

                    # inserisci downstream flanking se creato
                    if new_dw_id != 'NONE':
                        segs.insert(idx + 1, (new_dw_id, '+'))

                # Rimuovi i vecchi nodi dalle walk e dal grafo
                deleted_nodes = []
                for hap in selected_repetition[1]:
                    for n in hap[1]:
                        if n in nodes:
                            del nodes[n]
                            deleted_nodes.append(n)
                        # rimuovi ogni occorrenza in quella walk
                        wid = hap[0]
                        walks[wid]['segments'] = [(sid, o) for (sid, o) in walks[wid]['segments'] if sid != n]

                # Pulisci i vecchi link
                links = [lk for lk in links if lk[0] not in deleted_nodes and lk[2] not in deleted_nodes]

# -------------------- Scrittura GFA --------------------
with open(output_gfa, 'w') as gfa:
    # Header
    for header_line in header:
        gfa.write(f"{header_line}\n")

    # Segments
    for node_id, sequence in nodes.items():
        gfa.write(f"S\t{node_id}\t{sequence}\n")

    # Walks (W)
    for w_id, w in walks.items():
        walk_str = format_walk_segments(w['segments'])
        gfa.write(f"W\t{w['sample']}\t{w['hap']}\t{w['seqid']}\t{w['start']}\t{w['end']}\t{walk_str}\n")

    # Links
    for fr, fo, to, to_o, ov in links:
        gfa.write(f"L\t{fr}\t{fo}\t{to}\t{to_o}\t{ov}\n")

print(f"New GFA written to {output_gfa}")
