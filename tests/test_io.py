import gzip
from io import StringIO

import pytest

from grodon.io import highly_expressed_mask, read_coverage, read_fasta


def test_fasta_descriptions_case_and_selection():
    genes = read_fasta(
        StringIO(">a 30S ribosomal protein S1\natg\naaa\n>b ribosomal protein methylase\nccc\n")
    )
    assert tuple(genes) == ("ATGAAA", "CCC")
    assert genes.ids == ("a", "b")
    assert highly_expressed_mask(genes) == [True, False]
    assert highly_expressed_mask(genes, ids=["b"]) == [False, True]
    with pytest.raises(ValueError, match="missing"):
        highly_expressed_mask(genes, ids=["x"])


@pytest.mark.parametrize(
    "text", ["", "AAA\n", ">\nAAA\n", ">a\n", ">a\n>a\nAAA", ">a\nAAA\n>a\nTTT"]
)
def test_malformed_fasta(text):
    with pytest.raises(ValueError):
        read_fasta(StringIO(text))


def test_gzip_fasta_and_reordered_coverage(tmp_path):
    fasta = tmp_path / "genes.fa.gz"
    with gzip.open(fasta, "wt") as stream:
        stream.write(">a\nAAA\n>b\nCCC\n")
    genes = read_fasta(fasta)
    table = tmp_path / "depth.tsv"
    table.write_text("id\tcoverage\nb\t3\na\t1\n")
    assert read_coverage(table, genes) == [1.0, 3.0]


@pytest.mark.parametrize(
    "content",
    [
        "id\tcoverage\na\t1\n",
        "id\tcoverage\na\t1\na\t2\nb\t1\n",
        "id\tcoverage\na\tnan\nb\t1\n",
        "gene\tdepth\na\t1\n",
    ],
)
def test_bad_coverage(tmp_path, content):
    genes = read_fasta(StringIO(">a\nAAA\n>b\nCCC\n"))
    table = tmp_path / "coverage.tsv"
    table.write_text(content)
    with pytest.raises(ValueError):
        read_coverage(table, genes)
