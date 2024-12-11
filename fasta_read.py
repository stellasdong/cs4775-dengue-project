def parse_fasta(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as file:
        seq_id = None
        seq = ''
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    fasta_dict[seq_id] = seq
                seq_id = line[1:].split(' ')[0]
                seq = ''
            else:
                seq += line
                
        if seq_id:
            fasta_dict[seq_id] = seq
    
    return fasta_dict

fasta_file = 'test.fasta'
fasta_dict = parse_fasta(fasta_file)

# for seq_id, sequence in fasta_dict.items():
#     print(f"ID: {seq_id}")
#     print(f"Sequence: {sequence[:100]}...")

# print(fasta_dict)