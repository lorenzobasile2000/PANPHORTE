#!/usr/bin/env python3
"""
Script per raggruppare reads di un file FASTA per haplotype in batch di dimensione massima N (default 10), concatenandole e rinominandole.
Ogni nuovo read avrà header >R:<haplotype>:<numero_batch>
"""

import argparse

def parse_fasta(path):
    """
    Generatore che restituisce tuple (header, sequence) per ogni record FASTA
    """
    with open(path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header:
                    yield header, ''.join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield header, ''.join(seq_lines)

def main():
    parser = argparse.ArgumentParser(
        description='Raggruppa reads per haplotype in batch e concatena sequenze')
    parser.add_argument('-i', '--input', required=True,
                        help='Path del file FASTA in input')
    parser.add_argument('-o', '--output', required=True,
                        help='Path del file FASTA in output')
    parser.add_argument('-n', '--group-size', type=int, default=10,
                        help='Numero massimo di reads per batch (default: 10)')
    args = parser.parse_args()

    # Dizionario: chiave = haplotype (string), valore = lista di sequenze
    hap_reads = {}

    # Parsing del file FASTA
    for header, seq in parse_fasta(args.input):
        parts = header.split(':')
        # Ci aspettiamo formato R:<haplotype>:<read_index>
        if len(parts) >= 3 and parts[0] == 'R':
            hap = parts[1]
            hap_reads.setdefault(hap, []).append(seq)
        else:
            # Ignora record non conformi
            continue

    # Scrittura del FASTA di output
    with open(args.output, 'w') as out:
        # Ordina gli haplotype numericamente
        for hap in sorted(hap_reads.keys(), key=lambda x: int(x)):
            seq_list = hap_reads[hap]
            # Suddividi in batch di group_size
            for idx in range(0, len(seq_list), args.group_size):
                batch = seq_list[idx: idx + args.group_size]
                merged_seq = ''.join(batch)
                batch_num = idx // args.group_size + 1
                new_header = f">R:{hap}:{batch_num}"
                out.write(new_header + '\n')
                # Scrive la sequenza intera su una sola riga, senza wrapping
                out.write(merged_seq + '\n')

if __name__ == '__main__':
    main()
