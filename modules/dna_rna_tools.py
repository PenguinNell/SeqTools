def get_nucleic_acid_type(sequence,
                          alphabet_dna={'A', 'T', 'G', 'C'},
                          alphabet_rna={'A', 'U', 'G', 'C'}):

    unique_chars = set(sequence.upper())

    if unique_chars <= alphabet_dna:
        return 'DNA'
    elif unique_chars <= alphabet_rna:
        return 'RNA'
    else:
        return None


def is_nucleic_acid(sequence):

    seq_type = get_nucleic_acid_type(sequence)
    return seq_type in ['DNA', 'RNA']


def transcribe(sequence,
               mrna_compl={'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G',
                           'a': 'a', 't': 'u', 'c': 'c', 'g': 'g'}):

    seq_type = get_nucleic_acid_type(sequence)

    if seq_type == 'DNA':
        transcribe_seq = ''
        for nt in sequence:
            rna_nt = mrna_compl[nt]
            transcribe_seq += rna_nt
        return transcribe_seq

    elif seq_type == 'RNA':
        transcribe_seq = sequence
        return transcribe_seq
    else:
        print("Warning:", sequence, "is not a nucleotide sequence! Return None")
        return None


def reverse(sequence):

    if is_nucleic_acid(sequence):
        reverse_seq = sequence[::-1]
        return reverse_seq

    else:
        print("Warning:", sequence, "is not a nucleotide sequence! Return None")
        return None


def get_complement(sequence,
                   dna_compl={'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                              'a': 't', 't': 'a', 'g': 'c', 'c': 'g'},
                   rna_compl={'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                              'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}):

    seq_type = get_nucleic_acid_type(sequence)

    if not seq_type:
        print("Warning:", sequence, "is not a nucleotide sequence! Return None")
        return None

    else:
        if not (set(sequence.upper()) & {'T', 'U'}):
            print("Warning: sequence", sequence, "lacks T/U – processed as DNA by default")

        if seq_type == 'DNA':
            complement_seq = ''
            for nt in sequence:
                dna_nt = dna_compl[nt]
                complement_seq += dna_nt
            return complement_seq

        elif seq_type == 'RNA':
            complement_seq = ''
            for nt in sequence:
                rna_nt = rna_compl[nt]
                complement_seq += rna_nt
            return complement_seq


def get_reverse_complement(sequence):

    complement_seq = get_complement(sequence)
    if complement_seq is not None:
        reverse_complement_seq = reverse(complement_seq)
        return reverse_complement_seq
    else:
        return None


def run_dna_rna_tools(*args,
                      func_dict={'is_nucleic_acid': is_nucleic_acid,
                                 'transcribe': transcribe,
                                 'reverse': reverse,
                                 'complement': get_complement,
                                 'reverse_complement': get_reverse_complement}):

    if not args:
        print("Please provide arguments!")
        return None

    else:
        procedure = args[-1]
        seqs = args[0:-1]

        if not seqs and procedure in func_dict:
            print("Looks like you forgot to provide a nucleotide sequence(s). Please try again!")
            return None

        else:
            if procedure in func_dict:
                selected_function = func_dict[procedure]

                if len(seqs) == 1:
                    res = selected_function(seqs[0])
                    return res
                else:
                    res = []
                    for sq in seqs:
                        res.append(selected_function(sq))
                    return res

            else:
                print("Looks like you forgot to specify what to do with the sequence(s). Please try again!")
                return None
