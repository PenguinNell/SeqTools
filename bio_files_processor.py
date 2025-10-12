import os

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = "") -> None:
    """
    Converts a multi-line sequences in FASTA file to a single-line format for each sequence.

    Parameters
    ----------
    input_fasta : str
        path to the input multi-line FASTA file to be converted
    output_fasta : str, optional
        filename for the output single-line FASTA file.
        If not provided, a new file with '_oneline.fasta' will be created

    Returns
    -------
    None
        This function does not return any value. The processed sequences will be written to a FASTQ file.
    """

    if not os.path.isfile(input_fasta):
        print(f'File "{input_fasta}" does not exist or cannot be found. Please specify the full path')
        return None

    directory = os.path.dirname(input_fasta)

    if output_fasta == "":
        prefix = os.path.splitext(os.path.basename(input_fasta))[0]
        output_fasta = prefix + '_oneline.fasta'

    path_output_fasta = os.path.join(directory, os.path.basename(output_fasta))

    if os.path.isfile(path_output_fasta):
        open(path_output_fasta, 'w').close()

    with (
        open(input_fasta, 'r') as reads,
        open(path_output_fasta, 'a') as output_reads
    ):

        for line in reads:

            if line.startswith('>'):
                if 'read' not in locals():
                    read = [line.strip(), ""]
                else:
                    output_reads.write(read[0] + '\n')
                    output_reads.write(read[1] + '\n')
                    read = [line.strip(), ""]
            else:
                read[1] += line.strip()

    return None


def parse_blast_output(input_file: str, output_file: str = ""):
    """
     Extracts the top matched protein description for each query from a BLAST output file and saves the sorted list.

    Parameters
    ----------
    input_file : str
        path to the input BLAST results text file
    output_file : str, optional
        filename for the output file with the list of proteins.
        If not provided, a new file with '_proteins_list.txt' will be created

    Returns
    -------
    None
        The function writes the sorted list of top matched proteins to the output file.
    """

    if not os.path.isfile(input_file):
        print(f'File "{input_file}" does not exist or cannot be found. Please specify the full path')
        return None

    directory = os.path.dirname(input_file)

    if output_file == "":
        prefix = os.path.splitext(os.path.basename(input_file))[0]
        output_file = prefix + '_proteins_list.txt'

    path_output_file = os.path.join(directory, os.path.basename(output_file))

    proteins = []

    with open(input_file, 'r') as file:

        for line in file:

            if line.startswith('Sequences producing significant alignments:'):
                counter = 0

            if "counter" in locals():
                counter += 1

                if counter == 4:
                    proteins.append(line[:66].strip())

    proteins.sort(key=str.lower)

    with open(path_output_file, 'w') as output_file:
        output_file.write('\n'.join(proteins))

    return None


def select_genes_from_gbk_to_fasta(input_gbk,
                                   genes: str | list[str],
                                   n_before: int = 1,
                                   n_after: int = 1,
                                   output_fasta: str = ""):

    files = os.listdir(input_gbk)
    gbk_files = [f for f in files if f.endswith('.gbk')]

    for gbk_file in gbk_files:
        all_genes = []
        all_proteins = []

        with open(os.path.join(input_gbk, gbk_file), 'r') as file:

            protein_seq = ""
            gene_name = ""

            for line in file:
                if line.find('gene') != -1:
                    gene_name = line.strip()
                    all_genes.append(gene_name[gene_name.find('"'):].strip('"'))

                if line.find('translation') != -1:
                    protein_seq = ""

                if gene_name and line.find('CDS') != -1:
                    protein_seq = protein_seq[protein_seq.find('"'):].strip('"')
                    all_proteins.append(protein_seq)

                    protein_seq = ""

                protein_seq += line.strip()

        idx = all_genes.index(genes) # пока только со строкой

        neighbors = range(max(idx - n_before, 0), min(idx + n_after + 1, len(all_genes)))
        for i in neighbors:
            if i == idx:
                continue

            # print(f'>{all_genes[i]}\n{all_proteins[i]}')

# select_genes_from_gbk_to_fasta(input_gbk = ".", genes = "pxpB")
parse_blast_output("example_blast_results.txt")