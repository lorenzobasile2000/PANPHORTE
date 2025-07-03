import os

def process_fasta(file_path):
   
    if not os.path.exists(file_path):
        print(f"Error: file {file_path} doesn't exist.")
        return
    
    lengths = []
    
    with open(file_path, 'r') as file:
        read_name = None
        sequence = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                
                if read_name and sequence:
                    save_read(read_name, sequence)
                    lengths.append(f"{read_name}: {len(''.join(sequence))}\n")
                
            
                read_name = line[1:].split()[0]  
                sequence = []
            else:
                sequence.append(line.upper())
        
      
        if read_name and sequence:
            save_read(read_name, sequence)
            lengths.append(f"{read_name}: {len(''.join(sequence))}\n")
    
   
    with open("lengths.txt", 'w') as length_file:
        length_file.writelines(lengths)
    print("file created: lengths.txt")

def save_read(read_name, sequence):

    file_name = f"{read_name}.txt"
    with open(file_name, 'w') as out_file:
        out_file.write(''.join(sequence) + '\n')
    print(f"file created: {file_name}")


file_fasta = "chm13v2.0.fa"  
process_fasta(file_fasta)
