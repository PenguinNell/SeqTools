from __future__ import annotations # for classes that defined below

import os
from abc import ABC, abstractmethod
from typing import Self  # return the same type of object

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract class for biological sequences
    Provides interface for:
    - length calculation
    - indexing and slicing
    - print
    - alphabet check

    Parameters
    ----------
    seq : str
        Sequence string

    Raises
    ------
    ValueError
        If the sequence contains invalid characters
    """

    def __init__(self, seq: str) -> None:
        self.seq = str(seq)
        self.check_alphabet()

    def __len__(self) -> int:
        """
        Return sequence length

        Returns
        -------
        int
            Length of the stored sequence
        """
        return len(self.seq)

    def __getitem__(self, key: int | slice) -> str | Self:
        """
        Get element(s) by index or slice

        Parameters
        ----------
        key : int | slice
            Index (int) returns a single character (str)
            Slice returns an object of the same class as self

        Returns
        -------
        str | Self
            Single character for int index, or a new sequence string for slice
        """
        part = self.seq[key]
        if isinstance(key, slice):
            return self.__class__(part) # to keep data type (DNASequence -> DNASequence)
        return part

    def __str__(self) -> str:
        """
        Allows to print the sequence

        Returns
        -------
        str
            Stored sequence string
        """
        return self.seq

    @abstractmethod
    def check_alphabet(self) -> None:
        """
        Validate that sequence contains only allowed characters

        Raises
        ------
        ValueError
            If sequence contains invalid characters
        """
        raise NotImplementedError


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleotide sequences (DNA/RNA).

    Implements:
    - alphabet validation (based on subclass-defined alphabet)
    - complement / reverse / reverse_complement operations

    Notes
    -----
    Subclasses must define:
    - _ALPHABET : set[str]
    - _COMPLEMENT : dict[str, str]
    """

    _ALPHABET: set[str] | None = None
    _COMPLEMENT: dict[str, str] | None = None

    def check_alphabet(self) -> None:
        """
        Validate nucleotide alphabet

        Raises
        ------
        NotImplementedError
            For the instance of NucleicAcidSequence class
        ValueError
            If sequence contains invalid letters
        """
        if self._ALPHABET is None:
            raise NotImplementedError("Alphabet is not defined, use DNASequence or RNASequence class")

        invalid_letters = set(self.seq.upper()) - self._ALPHABET
        if invalid_letters:
            invalid_str = ", ".join(invalid_letters)
            raise ValueError(f"Invalid  letter(s) {invalid_str} in nucleotide sequence {self.seq}")

    def complement(self) -> Self:
        """
        Return the complement sequence

        Returns
        -------
        Self
            Complement sequence of the same class as self

        Raises
        ------
        NotImplementedError
            For the instance of NucleicAcidSequence class
        """
        if self._COMPLEMENT is None:
            raise NotImplementedError("Complement is not defined, use DNASequence or RNASequence class")

        output = []
        for nt in self.seq:
            output.append(self._COMPLEMENT[nt])
        return self.__class__("".join(output)) # keep class type

    def reverse(self) -> Self:
        """
        Return reversed sequence

        Returns
        -------
        Self
            Reversed sequence of the same class as self
        """
        return self.__class__(self.seq[::-1])

    def reverse_complement(self) -> Self:
        """
        Return reverse complement sequence

        Returns
        -------
        Self
            Reverse complement sequence of the same class as self
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    DNA sequence class
    """

    _ALPHABET = set("ATGC")
    _COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                   'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    def transcribe(self) -> RNASequence:
        """
        Transcribe DNA to RNA (T -> U), keep character case

        Returns
        -------
        RNASequence
            RNA sequence produced by transcription
        """
        output = []
        for nt in self.seq:
            if nt == "T":
                output.append("U")
            elif nt == "t":
                output.append("u")
            else:
                output.append(nt)
        return RNASequence("".join(output))


class RNASequence(NucleicAcidSequence):
    """
    RNA sequence class
    """

    _ALPHABET = set("AUGC")
    _COMPLEMENT = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                   'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


class AminoAcidSequence(BiologicalSequence):
    """
    Amino acid sequence class
    """

    _ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self) -> None:
        """
        Validate amino acid alphabet

        Raises
        ------
        ValueError
            If sequence contains invalid letters
        """
        invalid_letters = set(self.seq.upper()) - self._ALPHABET
        if invalid_letters:
            invalid_str = ", ".join(invalid_letters)
            raise ValueError(f"Invalid letter(s) {invalid_str} in amino acid sequence {self.seq}")

    def aa_frequencies(self) -> dict[str, float]:
        """
        Compute relative frequencies of amino acids in the sequence

        Returns
        -------
        dict[str, float]
            Dictionary for mapping letters and their frequency, rounded to 2 decimals
        """
        counts = {}
        for aa in self.seq.upper():
            if aa in counts:
                counts[aa] += 1
            else:
                counts[aa] = 1

        freqs = {}
        for aa, count in counts.items():
            freqs[aa] = round(count / len(self.seq), 2)

        return freqs


def filter_fastq(input_fastq: str,
                 output_fastq: str = "",
                 overwrite: bool = False,
                 gc_bounds: int | float | tuple = (0, 100),
                 length_bounds: int | float | tuple = (0, 2**32),
                 quality_threshold: int | float = 0) -> None:
    """
    Filters FASTQ file by GC content, length, and quality

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
        This function does not return any value. The filtered data is written to a FASTQ file
    """

    if not os.path.isfile(input_fastq):
        print(f'File "{input_fastq}" does not exist or cannot be found. Please specify the full path')
        return None

    directory = os.path.dirname(input_fastq)

    filtered_folder = os.path.join(directory, 'filtered')
    if not os.path.exists(filtered_folder):
        os.makedirs(filtered_folder)

    if output_fastq == "":
        prefix = os.path.splitext(os.path.basename(input_fastq))[0]
        output_fastq = prefix + '_filtered.fastq'

    path_output_fastq = os.path.join(filtered_folder, os.path.basename(output_fastq))

    if os.path.isfile(path_output_fastq) and not overwrite:
        print(f'File "{output_fastq}" is exist! Use overwrite=True to overwrite it')
        return None

    with (
        open(input_fastq, 'r') as reads,
        open(path_output_fastq, 'w') as output_file
    ):
        for record in SeqIO.parse(reads, "fastq"):

            gc_percent = gc_fraction(record.seq) * 100
            if isinstance(gc_bounds, (int, float)):
                gc_bounds = (0, gc_bounds)
            if not (gc_bounds[0] <= gc_percent <= gc_bounds[1]):
                continue

            if isinstance(length_bounds, (int, float)):
                length_bounds = (0, length_bounds)
            if not (length_bounds[0] <= len(record) <= length_bounds[1]):
                continue

            phred = record.letter_annotations["phred_quality"]
            mean_q = sum(phred) / len(phred)
            if not (mean_q >= quality_threshold):
                continue

            SeqIO.write(record, output_file, "fastq")

    if os.path.getsize(path_output_fastq) == 0:
        print("No sequences passed the filters! Output file is empty")

    return None