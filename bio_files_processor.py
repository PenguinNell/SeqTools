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

convert_multiline_fasta_to_oneline(input_fasta = "example_multiline_fasta.fasta", output_fasta="res.fasta")
