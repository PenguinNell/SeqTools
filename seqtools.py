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

    Arguments:
    *args: any number of string arguments, where:
        - the first N arguments are nucleotide sequences
        - the last argument is the procedure to apply
          (supported procedures: 'is_nucleic_acid', 'transcribe', 'reverse',
          'complement', 'reverse_complement')

    Returns result of the specified procedure: str / list
    If the procedure is not supported - returns None
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


def filter_fastq(fastq_seqs: dict,
                 gc_bounds: int | float | tuple = (0, 100),
                 length_bounds: int | float | tuple = (0, 2**32),
                 quality_threshold: int | float = 0) -> dict | None:
    """
    Filters FASTQ sequences by GC content, length, and quality.

    Arguments:
    fastq_seqs (dict): dictionary with sequence names (IDs) as keys and tuples of (sequence, quality) as values
    gc_bounds (int | float | tuple, optional): GC content range in percent (default: (0, 100))
        Can be a tuple (min, max) or a single number (upper threshold)
    length_bounds (int | float | tuple, optional): length range (default: (0, 2**32))
         Can be a tuple (min, max) or a single number (upper threshold)
    quality_threshold (int | float, optional): minimum average quality (default: 0)

    Returns: dict | None (filtered sequences, or None if none pass)
    """

    filtered_fastq_seqs = {}

    for seq_id, (seq, quality) in fastq_seqs.items():
        if (
                fastq_filter.is_in_gc_bounds(seq, gc_bounds) and
                fastq_filter.is_in_length_bounds(seq, length_bounds) and
                fastq_filter.is_above_quality_threshold(quality, quality_threshold)
        ):
            filtered_fastq_seqs[seq_id] = (seq, quality)

    if not filtered_fastq_seqs:
        print("No sequences passed the filters! Return None")
        return None

    return filtered_fastq_seqs