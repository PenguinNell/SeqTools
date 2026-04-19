import pytest

from seqtools import DNASequence, RNASequence, AminoAcidSequence, filter_fastq
from bio_files_processor import convert_multiline_fasta_to_oneline


class TestDNASequence:
    def test_getitem_index_return(self):
        seq = DNASequence("AGATCTGGAATGAAC")
        assert isinstance(seq[0], str)
        assert seq[0] == "A"
        assert seq[-1] == "C"

    def test_getitem_slice_return(self):
        seq = DNASequence("AGATCTGGAATGAAC")
        assert isinstance(seq[3:7], DNASequence)
        assert str(seq[3:7]) == "TCTG"
        assert str(seq[::-1]) == "CAAGTAAGGTCTAGA"

    def test_reverse_complement(self):
        seq = DNASequence("AGATCTGGAATGAAC")
        assert str(seq.reverse_complement()) == "GTTCATTCCAGATCT"

    def test_transcribe_case_match(self):
        seq = DNASequence("AgATCtggaATGaAC")
        rna = seq.transcribe()
        assert isinstance(rna, RNASequence)
        assert str(rna) == "AgAUCuggaAUGaAC"

    def test_invalid_seq(self):
        with pytest.raises(ValueError):
            DNASequence("AGADWVGЯЖATGAAC")


class TestAminoAcidSequence:
    def test_aa_frequencies(self):
        seq = AminoAcidSequence("MHlFirlQfWTkFVaALCnY")
        expected_result = {"A": 0.1, "C": 0.05, "F": 0.15, "H": 0.05, "I": 0.05, "K": 0.05, "L": 0.15, "M": 0.05,
                           "N": 0.05, "Q": 0.05, "R": 0.05, "T": 0.05, "V": 0.05, "W": 0.05, "Y": 0.05}
        assert seq.aa_frequencies() == pytest.approx(expected_result, abs=1e-12)


class TestBioFileProcessing:
    def test_filter_fastq_creates_filtered_file_and_keeps_passing_read(self, tmp_path):
        inp = tmp_path / "reads.fastq"
        inp.write_text(
            "@read1\n"
            "TAATTTATATGACAG\n" # GC 20%
            "+\n"
            "IIIIIIIIIIIIIII\n" # 40
            "@read2\n"
            "TTAATCTC\n" # GC 25%
            "+\n"
            "IIHGDCA@\n" # 36.13
            "@read3\n"
            "GTTGGTGTAGGTCCCGACCGACGCG\n" # CG 68%
            "+\n"
            "IIIIIIIIHGFEDCBA@?>=<;:98\n" # 35.12
            "@read4\n"
            "CGCCTCGTTGG\n" # GC 72.7%
            "+\n"
            "55555555555\n" # 20
        )

        filter_fastq(
            input_fastq=str(inp),
            output_fastq="out1.fastq",
            overwrite=True,
            gc_bounds=50,
            quality_threshold=38
        )

        out1_path = tmp_path / "filtered" / "out1.fastq"
        assert out1_path.exists()

        out1 = out1_path.read_text()
        assert "@read1" in out1
        assert all(read not in out1 for read in ["@read2", "@read3", "@read4"])

        filter_fastq(
            input_fastq=str(inp),
            output_fastq="out2.fastq",
            overwrite=True,
            gc_bounds=(50, 100),
            quality_threshold=30
        )

        out2_path = tmp_path / "filtered" / "out2.fastq"
        assert out2_path.exists()

        out2 = out2_path.read_text()
        assert "@read3" in out2
        assert all(read not in out2 for read in ["@read1", "@read2", "@read4"])

    def test_convert_multiline_fasta_to_oneline(self, tmp_path):
        inp = tmp_path / "input.fasta"
        out = tmp_path / "output.fasta"

        inp.write_text(
            ">seq1\n"
            "AGATC\n"
            "TGGAAT\n"
            "GAAC\n"
            ">seq2\n"
            "A\n"
            "A\n"
            "T\n"
        )

        convert_multiline_fasta_to_oneline(str(inp), str(out))

        assert out.exists()
        assert out.read_text() == (
            ">seq1\n"
            "AGATCTGGAATGAAC\n"
            ">seq2\n"
            "AAT\n"
        )