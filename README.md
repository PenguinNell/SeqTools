# SeqTools
**SeqTools** is a Python toolkit designed for working with nucleotide sequences. It provides essential functionality for DNA/RNA sequence operations, FASTQ filtering, and BLAST and FASTA files processing.


## Features
- **DNA/RNA Tools**: Validate nucleotide sequences or perform operations like transcription, reverse, complement, and reverse_complement.
- **Sequence Filtering**: Filter FASTQ files by GC content, length, or base quality. Filtered file will be saved in 'filtered' subdirectory. 
  - Filtering parameters are **inclusive** (e.g. `gc_bounds=(40, 60)` means that sequences with a GC content of **40% and 60%** will also be included). 
  - If only one number is provided for `gc_bounds` or `length_bounds`, it is treated as the **upper threshold**, with min defaulting to 0.
- **Multi-line FASTA Conversion**: Convert multi-line FASTA sequences to single-line format for easier downstream processing.
- **BLAST Parsing**: Extract the top-matched protein description for each query from BLAST text output and save a sorted protein list.


## Installation
To use SeqTools, clone the repository:
```bash
git clone https://github.com/PenguinNell/SeqTools.git
```

## Usage

### Example 1: Using DNA/RNA tools

```python
from seqtools import run_dna_rna_tools

result = run_dna_rna_tools('ATG', 'transcribe')
print(result)  # Output: 'AUG' (RNA sequence)
```

### Example 2: Filtering FASTQ files

```python
from seqtools import filter_fastq

filter_fastq(
    input_fastq="file.fastq",
    output_fastq="output_filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(5, 10),
    quality_threshold=30
)
# Output: Filtered FASTQ file in 'filtered' subdirectory
```

### Example 3: Converting Multi-line FASTA to One-line Format
```python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline('input.fasta', 'output_oneline.fasta')
# Output: Each sequence in 'output_oneline.fasta' will be in single-line format
```

### Example 4: Parsing BLAST Output for Top Protein Hits
```python
from bio_files_processor import parse_blast_output

parse_blast_output('results.txt', 'proteins.txt')
# Output: Each line in 'proteins.txt' contains the best protein match for each query, sorted alphabetically
```