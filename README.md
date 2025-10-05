# SeqTools
**SeqTools** is a Python toolkit designed for working with nucleotide sequences. It provides essential functionality for DNA/RNA sequence operations and FASTQ filtering.


## Features
- **DNA/RNA Tools**: Validate nucleotide sequences or perform operations like transcription, reverse, complement, and reverse_complement.
- **Sequence Filtering**: Filter FASTQ sequences by GC content, length, or base quality. 
  - Filtering parameters are **inclusive** (e.g. `gc_bounds=(40, 60)` means that sequences with a GC content of **40% and 60%** will also be included). 
  - If only one number is provided for `gc_bounds` or `length_bounds`, it is treated as the **upper threshold**, with min defaulting to 0.


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

### Example 2: Filtering FASTQ sequences

```python
from seqtools import filter_fastq

# Input FASTQ data
fastq_seqs = {'@seq1': ("ATGCATGC", "IIIIIIII")}

# Apply filters
filtered = filter_fastq(
    fastq_seqs=fastq_seqs,
    gc_bounds=(40, 60),
    length_bounds=(5, 10),
    quality_threshold=30
)
print(filtered)  # Output: Filtered sequence dictionary
```