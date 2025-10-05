import modules.dna_rna_tools as nuclei
import modules.fastq_filter_checkers as fastq_filter

func_dict={'is_nucleic_acid': nuclei.is_nucleic_acid,
           'transcribe': nuclei.transcribe,
           'reverse': nuclei.reverse,
           'complement': nuclei.get_complement,
           'reverse_complement': nuclei.get_reverse_complement}


def run_dna_rna_tools(*args):
    *seqs, procedure = args

    if not seqs and procedure in func_dict:
        print("Looks like you forgot to provide a nucleotide sequence(s). Please try again!")
        return None

    if procedure in func_dict:
        selected_function = func_dict[procedure]

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
                 quality_threshold: int | float = 0) -> dict | None: # по заданию первый аргумент должен называться seqs, но я решила дать ему другое название, чтобы не было путанницы со списком сиквенсов НЕ в fastq-формате
    # check no input
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