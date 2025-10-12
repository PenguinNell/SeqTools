import os
import modules.dna_rna_tools as nuclei
import modules.fastq_filter_checkers as fastq_filter


dna_rna_func={'is_nucleic_acid': nuclei.is_nucleic_acid,
              'transcribe': nuclei.transcribe,
              'reverse': nuclei.reverse,
              'complement': nuclei.get_complement,
              'reverse_complement': nuclei.get_reverse_complement}


def run_dna_rna_tools(*args: str) -> str | list | None:
    """
    Applies a specified procedure to one or multiple nucleotide sequences

    Parameters
    ----------
    *args :
        any number of string arguments, where: the first N arguments are nucleotide sequences
        and the last argument is the procedure to apply
        (supported procedures: 'is_nucleic_acid', 'transcribe', 'reverse',
        'complement', 'reverse_complement')

    Returns
    -------
    str | list | None
        Returns result of the specified procedure (str or list). If the procedure is not supported - returns None
    """

    *seqs, operation = args

    if not seqs and operation in dna_rna_func:
        print("Looks like you forgot to provide a nucleotide sequence(s). Please try again!")
        return None

    if operation in dna_rna_func:
        selected_function = dna_rna_func[operation]

        if len(seqs) == 1:
            res = selected_function(seqs[0])
            return res

        res = []
        for sq in seqs:
            res.append(selected_function(sq))
        return res

    print("Looks like you forgot to specify what to do with the sequence(s). Please try again!")
    return None


def filter_fastq(input_fastq: str,
                 output_fastq: str = None,
                 overwrite: bool = False,
                 gc_bounds: int | float | tuple = (0, 100),
                 length_bounds: int | float | tuple = (0, 2**32),
                 quality_threshold: int | float = 0) -> None:
    """
    Filters FASTQ file by GC content, length, and quality.

    Parameters
    ----------
    input_fastq : str
        path to the input FASTQ file to be filtered
    output_fastq : str, optional
        name to the output filtered FASTQ file, will be created in a 'filtered' subdirectory
    overwrite : bool, optional
        if True, overwrites the output file if it exists
    gc_bounds : int | float | tuple, optional
        GC content range in percent (default: (0, 100)). Can be a tuple (min, max) or a single number (upper threshold)
    length_bounds  : int | float | tuple, optional
        length range (default: (0, 2**32)). Can be a tuple (min, max) or a single number (upper threshold)
    quality_threshold : int | float, optional
        minimum average quality (default: 0)

    Returns
    -------
    None
        This function does not return any value. The filtered data is written to a FASTQ file.
    """

    if not os.path.isfile(input_fastq):
        print(f'File "{input_fastq}" does not exist or cannot be found. Please specify the full path')
        return None

    directory = os.path.dirname(input_fastq)

    if output_fastq == "":
        prefix = os.path.splitext(os.path.basename(input_fastq))[0]
        output_fastq = prefix + '_filtered.fastq'

    path_output_fastq = os.path.join(directory, 'filtered', os.path.basename(output_fastq))

    if path_output_fastq:
        print(f'File "{output_fastq}" is exist! Use overwrite=True to overwrite it')
        return None

    if overwrite:
        open(path_output_fastq, 'w').close()

    with (
        open(input_fastq, 'r') as reads,
        open(path_output_fastq, 'a') as output_file
    ):
        read = []

        for line in reads:
            read.append(line.strip())

            if len(read) == 4:
                seq_id = read[0]
                seq = read[1]
                separator = read[2]
                quality = read[3]

                if fastq_filter.is_read_good(seq, quality, gc_bounds, length_bounds, quality_threshold):
                    output_file.write(seq_id + '\n')
                    output_file.write(seq + '\n')
                    output_file.write(separator + '\n')
                    output_file.write(quality + '\n')

                read = []

    return None