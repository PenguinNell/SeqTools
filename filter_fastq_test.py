import pytest
from pathlib import Path
from Bio import SeqIO

from seqtools import filter_fastq


@pytest.fixture
def test_fastq_file(tmp_path):
    inp = tmp_path / "input.fastq"
    inp.write_text(
        "@read1\n"
        "TAATTTATATGACAG\n" # len 15, GC 20%
        "+\n"
        "IIIIIIIIIIIIIII\n" # 40
        "@read2\n"
        "TTAATCTC\n" # len 8, GC 25%
        "+\n"
        "IIHGDCA@\n" # 36.13
        "@read3\n"
        "GTTGGTGTAGGTCCCGACCGACGCG\n" # len 25, CG 68%
        "+\n"
        "IIIIIIIIHGFEDCBA@?>=<;:98\n" # 33.88
        "@read4\n"
        "CGCCTCGTTGG\n" # len 11, GC 72.7%
        "+\n"
        "55555555555\n", # 20
        encoding="utf-8"
    )
    return inp


def _create_output_path(input_fastq, output_fastq=""):
    inp = Path(input_fastq)
    out_dir = inp.parent / "filtered"
    prefix = inp.stem
    if output_fastq == "":
        out_file = out_dir / f"{prefix}_filtered.fastq"
    else:
        out_file = out_dir / output_fastq
    return out_file


class TestCheckInput:
    def test_no_input_file(self, tmp_path, capsys):
        pseudo_file = tmp_path / "this_file_does_not_exist.fastq"
        filter_fastq(pseudo_file)
        captured = capsys.readouterr()
        assert f'File "{pseudo_file}" does not exist or cannot be found. Please specify the full path' in captured.out

    def test_invalid_fastq_error(self, tmp_path):
        bad_file = tmp_path / "bad.fastq"
        bad_file.write_text(
            "@bad_read1\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIII\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError):
            filter_fastq(bad_file)


class TestCheckParams:
    def test_gc_bounds_upper(self, test_fastq_file):
        filter_fastq(test_fastq_file, gc_bounds=25)
        out_file = _create_output_path(test_fastq_file)
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read1", "read2"]

    def test_gc_bounds_tuple(self, test_fastq_file):
        filter_fastq(test_fastq_file, gc_bounds=(25, 70))
        out_file = _create_output_path(test_fastq_file)
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read2", "read3"]

    def test_len_bounds_tuple(self, test_fastq_file):
        filter_fastq(test_fastq_file, length_bounds=(10, 20))
        out_file = _create_output_path(test_fastq_file)
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read1", "read4"]

    def test_quality_threshold_filter(self, test_fastq_file):
        filter_fastq(test_fastq_file, quality_threshold=35)
        out_file = _create_output_path(test_fastq_file)
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read1", "read2"]


class TestCheckOutput:
    def test_creates_filtered_dir_and_def_output(self, test_fastq_file):
        filter_fastq(test_fastq_file)  # all should pass
        out_file = _create_output_path(test_fastq_file)
        assert out_file.exists()
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read1", "read2", "read3", "read4"]

    def test_creates_given_output(self, test_fastq_file):
        filter_fastq(input_fastq=test_fastq_file, output_fastq="test.fastq")  # all should pass
        out_file = _create_output_path(input_fastq=test_fastq_file, output_fastq="test.fastq")
        assert out_file.exists()
        ids = [read.id for read in SeqIO.parse(out_file, "fastq")]
        assert ids == ["read1", "read2", "read3", "read4"]

    def test_not_overwrite_output(self, test_fastq_file, capsys):
        # 1st call
        filter_fastq(input_fastq=test_fastq_file, output_fastq="test.fastq")  # all should pass
        out_file1 = _create_output_path(input_fastq=test_fastq_file, output_fastq="test.fastq")
        ids = [read.id for read in SeqIO.parse(out_file1, "fastq")]

        # 2nd call
        filter_fastq(input_fastq=test_fastq_file,
                     output_fastq="test.fastq",
                     length_bounds=(10, 20))  # only @read1 and @read4 should pass
        captured = capsys.readouterr()
        assert 'Use overwrite=True to overwrite it' in captured.out
        out_file2 = _create_output_path(input_fastq=test_fastq_file, output_fastq="test.fastq")
        assert ids == [read.id for read in SeqIO.parse(out_file2, "fastq")]

    def test_empty_output(self, test_fastq_file, capsys):
        filter_fastq(test_fastq_file, gc_bounds=0)
        out_file = _create_output_path(test_fastq_file)
        assert out_file.is_file()
        assert out_file.stat().st_size == 0

        captured = capsys.readouterr()
        assert "No sequences passed the filters! Output file is empty" in captured.out