import re
import subprocess
import json

# Parameters
gfa_file = "example3.gfa"
out_json = "OUT_JSON"
fasta = "FASTA"
file_path = "OUT_JSON"


def read_gfa(file_gfa):
    nodes = {}

    with open(file_gfa, 'r') as gfa:
        for line in gfa:
            line = line.strip()
            if not line:
                continue

            tipo_record = line[0]

            if tipo_record == 'S':
                seq = line.split('\t')
                if len(seq) >= 3:
                    id_node = seq[1]
                    sequence = seq[2]
                    nodes[id_node] = sequence  # Aggiungi il nodo e la sua sequenza

    return nodes


def regex(sequence):
    pattern = re.compile(r'(.+?)\1+')
    matches = pattern.finditer(sequence)

    tandem_repetitions = []
    for match in matches:
        repetition = match.group(1)
        times = len(match.group(0)) // len(repetition)
        position = match.start()
        tandem_repetitions.append((repetition, times, position))

    if tandem_repetitions:
        for repetition, times, position in tandem_repetitions:
            print(f"{repetition} repeated {times} times at position {position}\n")
    else:
        print("No tandem repetition find.\n")

    return tandem_repetitions


nodes = read_gfa(gfa_file)
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
        if len(bubble['inside']) > 1:

            for path in bubble['inside']:
                sequence = ''
                if isinstance(path, str):
                    sequence = sequence + nodes[path]
                else:
                    for node in path:
                        sequence = sequence + nodes[node]
                print(f"Calling RegexPy with: {sequence}")
                rep = regex(sequence)
        else:
            print("No further analysis needed!")
        print("-----")
