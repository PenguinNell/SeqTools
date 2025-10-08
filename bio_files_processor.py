def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = ""):

    if output_fasta == "":
        prefix = input_fasta.split('.')[0]
        output_fasta = prefix + '_oneline.fasta'

    oneline_reads = {}

    with open(input_fasta, 'r') as reads:
        read = [reads.readline().strip()]

        for line in reads:

            if line.startswith('>'):
                oneline_reads[read[0]] = ''.join(read[1:])
                read = [line.strip()]
            else:
                read.append(line.strip())

        oneline_reads[read[0]] = ''.join(read[1:])

    with open(output_fasta, 'w') as output_reads:
        for name, seq in oneline_reads.items():
            output_reads.write(name + '\n')
            output_reads.write(seq + '\n')


def parse_blast_output(input_file: str, output_file: str = ""):

    with open(input_file, 'r') as file:
        lines = file.readlines()

    first_lines = []

    for l in range(len(lines)):
        line = lines[l]

        if line.startswith('Sequences producing significant alignments:'):
            first_lines.append(lines[l+3].strip())

    proteins = []
    for l in range(len(first_lines)):
        first_line = first_lines[l]

        idx_spaces = first_line.find('   ')
        idx_dots = first_line.find('...')

        indices = [i for i in [idx_spaces, idx_dots] if i != -1]

        stop = min(indices)
        if stop == idx_dots:
            protein_name = first_line[:(stop+3)]
        else:
            protein_name = first_line[:stop]

        proteins.append(protein_name)

    proteins.sort(key=str.lower)

    if output_file == "":
        prefix = input_file.split('.')[0]
        output_file = prefix + '_proteins.txt'

    with open(output_file, 'w') as output_reads:
        for i in range(len(proteins)):
            output_reads.write(proteins[i] + '\n')