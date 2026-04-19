# SeqTools
**SeqTools** is a Python toolkit designed for working with biological sequences (DNA, RNA, protein). It provides essential functionality for DNA/RNA/amino acid sequence operations, FASTQ filtering, and BLAST and FASTA files processing.

## Features

### 1) Biological sequence classes (OOP)
Implemented in `seqtools.py`:

- **Common interface** (class `BiologicalSequence`)
  - `len(seq_obj)` support
  - indexing and slicing via `seq_obj[i]` and `seq_obj[i:j]`
  - printing via `print(seq_obj)`
  - alphabet validation on initialization (raises `ValueError`)

- **Nucleotide sequences** (classes `DNASequence` and `RNASequence`)
  - `.complement()`
  - `.reverse()`
  - `.reverse_complement()`
  - `DNASequence.transcribe()` returns an `RNASequence`
  - alphabet validation:
    - DNA: `A,T,G,C` (case-insensitive)
    - RNA: `A,U,G,C` (case-insensitive)

- **Protein sequences** (class `AminoAcidSequence`)
  - alphabet validation: 20 standard amino acids `ACDEFGHIKLMNPQRSTVWY` (case-insensitive)
  - `.aa_frequencies()` returns relative amino acid frequencies

### 2) FASTQ filtering (Biopython-based)
Implemented in `seqtools.py` as `filter_fastq`.

Filter FASTQ files by GC content, length, or base quality. Filtered file will be saved in 'filtered' subdirectory. 

Filter conditions:
- **GC content** in percent
- **sequence length**
- **mean Phred quality**

Notes:
- All bounds are **inclusive** (e.g. `gc_bounds=(40, 60)` means that sequences with a GC content of **40% and 60%** will also be included).
- If only one number is provided for `gc_bounds` or `length_bounds`, it is treated as the **upper threshold**, with min defaulting to 0.
- If no reads pass filters, the output file is created but empty and a message is printed.

### 3) FASTA and BLAST utilities
Implemented in `bio_files_processor.py`:

- **Multi-line FASTA Conversion**: Convert multi-line FASTA sequences to single-line format for easier downstream processing.
- **BLAST Parsing**: Extract the top-matched protein description for each query from BLAST text output and save a sorted protein list.

## Installation
To use SeqTools, clone the repository:
```bash
git clone https://github.com/PenguinNell/SeqTools.git
```

Install dependencies (Biopython is required for FASTQ filtering):
```
pip install -r requirements.txt
```

## Usage

### Example 1: Using biological sequence operations

```python
from seqtools import DNASequence, RNASequence, AminoAcidSequence

dna = DNASequence("AtGc")
print(dna)  # Output: AtGc
print(dna.complement())  # Output: TaCg
print(dna.reverse())  # Output: cGtA
print(dna.reverse_complement())  # Output: gCaT

rna = dna.transcribe()
print(rna)  # Output: AuGc
print(type(rna)) # Output: <class 'seqtools.RNASequence'>

protein = AminoAcidSequence("aAaCC")
print(protein.aa_frequencies())  # Output: {'A': 0.6, 'C': 0.4} 
```

### Example 2a : Filtering FASTQ files (Python API)

```python
from seqtools import filter_fastq

filter_fastq(
    input_fastq="file.fastq",
    output_fastq="output_filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(5, 10),
    quality_threshold=30
)
#  Output: Filtered FASTQ file in 'filtered' subdirectory
```

### Example 2b: Filtering FASTQ files (CLI)
A command-line interface for FASTQ filtering is available via the script `filter_fastq_cli.py`.

Show help:
```bash
python filter_fastq_cli.py --help
```

Minimal run (only input file; all filters are defaults):
```bash
python filter_fastq_cli.py reads.fastq
```

CLI supports the same filters:
```bash
python filter_fastq_cli.py reads.fastq --gc-bounds 40 60 --length-bounds 5 10 -q 30 --overwrite
```

If you want to provide only an upper bound, use the `--gc-upper` / `--length-upper` options instead of `--gc-bounds` / `--length-bounds`:
```bash
python filter_fastq_cli.py reads.fastq --gc-upper 60 --length-upper 10 -q 30
```
Use either `--gc-upper` or `--gc-bounds` (same for length), not both. 

### Example 3: Converting Multi-line FASTA to One-line Format

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline('input.fasta', 'output_oneline.fasta')
#  Output: Each sequence in 'output_oneline.fasta' will be in single-line format
```

### Example 4: Parsing BLAST Output for Top Protein Hits

```python
from bio_files_processor import parse_blast_output

parse_blast_output('results.txt', 'proteins.txt')
#  Output: Each line in 'proteins.txt' contains the best protein match for each query, sorted alphabetically
```

## Contact
For questions, suggestions, or contributions, please contact me at:
ansycheva26@gmail.com