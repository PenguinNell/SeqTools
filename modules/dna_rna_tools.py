alphabet_dna = {'A', 'T', 'G', 'C'}
alphabet_rna = {'A', 'U', 'G', 'C'}
dna_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
rna_complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                  'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
mrna_complement = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G',
                   'a': 'a', 't': 'u', 'c': 'c', 'g': 'g'}


def is_dna(sequence: str) -> bool:
    unique_chars = set(sequence.upper())

    if unique_chars <= alphabet_dna:
        return True

    return False


def  is_rna(sequence: str) -> bool:
    unique_chars = set(sequence.upper())

    if unique_chars <= alphabet_rna:
        return True

    return False


def is_nucleic_acid(sequence: str) -> bool:
    return is_dna(sequence) or is_rna(sequence)


def transcribe(sequence: str) -> str | None:
    if is_dna(sequence):
        transcribe_seq = ''.join([mrna_complement[nt] for nt in sequence])
        return transcribe_seq

    if is_rna(sequence):
        return sequence

    print("Warning:", sequence, "is not a nucleotide sequence! Return None")
    return None


def reverse(sequence: str) -> str | None:
    if is_nucleic_acid(sequence):
        return sequence[::-1]

    print("Warning:", sequence, "is not a nucleotide sequence! Return None")
    return None


def get_complement(sequence: str) -> str | None:
    if not is_nucleic_acid(sequence):
        print("Warning:", sequence, "is not a nucleotide sequence! Return None")
        return None

    if not (set(sequence.upper()) & {'T', 'U'}):
        print("Warning: sequence", sequence, "lacks T/U – processed as DNA by default")

    if is_dna(sequence):
        return ''.join([dna_complement[nt] for nt in sequence])

    if is_rna(sequence):
        return ''.join([rna_complement[nt] for nt in sequence])


def get_reverse_complement(sequence: str) -> str | None:
    complement_seq = get_complement(sequence)

    if complement_seq is not None:
        return reverse(complement_seq)

    return None