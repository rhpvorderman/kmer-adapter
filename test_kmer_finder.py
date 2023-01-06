import operator
import string

import pytest

from cutadapt._align import _upper_table, _acgt_table, _iupac_table
from kmer_finder import KmerFinder

UPPER_TABLE: bytes = _upper_table()
ACGT_TABLE: bytes = _acgt_table()
IUPAC_TABLE: bytes = _iupac_table()

KMER_FINDER_TESTS = [
    # kmer, start, stop, ref_wildcards, query_wildcards, sequence, expected
    ("ACGT", 0, None, False, False, "ACGTACG", True),
    ("ACGT", 0, None, False, False, "ACgtACG", True),
    ("acgt", 0, None, False, False, "ACgtACG", True),
    ("ACGT", 0, None, False, False, "acgtacg", True),
    ("ACGT", 0, None, False, False, "gacgact", False),
    ("ACGT", 0, None, False, True, "ACGNACG", True),
    ("ACGT", 0, None, False, False, "ACGNACG", False),
    ("ACGN", 0, None, True, False, "ACGTACG", True),
    ("ACGN", 0, None, True, False, "ACGxACG", True),
    ("ACKN", 0, None, True, False, "ACGTACG", True),
    ("ACKN", 0, None, True, True, "ACWRACG", True),
    ("ACKN", 0, None, True, True, "ACWxACG", False),
]


@pytest.mark.parametrize(
    [
        "kmer",
        "start",
        "stop",
        "ref_wildcards",
        "query_wildcards",
        "sequence",
        "expected",
    ],
    KMER_FINDER_TESTS,
)
def test_kmer_finder(
    kmer: str,
    start: int,
    stop: int,
    ref_wildcards: bool,
    query_wildcards: bool,
    sequence: str,
    expected: bool,
):
    kmer_finder = KmerFinder([(start, stop, [kmer])], ref_wildcards, query_wildcards)
    assert kmer_finder.kmers_present(sequence) is expected


@pytest.mark.parametrize(
    ["ref_table", "query_table", "comp_op", "ref_wildcards", "query_wildcards"],
    [
        (UPPER_TABLE, UPPER_TABLE, operator.eq, False, False),
        (IUPAC_TABLE, ACGT_TABLE, operator.and_, True, False),
        (ACGT_TABLE, IUPAC_TABLE, operator.and_, False, True),
        (IUPAC_TABLE, IUPAC_TABLE, operator.and_, True, True),
    ],
)
def test_kmer_finder_per_char_matching(
    ref_table, query_table, comp_op, ref_wildcards, query_wildcards
):
    for char in string.ascii_letters:
        kmer_finder = KmerFinder(
            [(0, None, [char])],
            ref_wildcards=ref_wildcards,
            query_wildcards=query_wildcards,
        )
        ref_char = ref_table[ord(char)]
        for comp_char in string.ascii_letters:
            query_char = query_table[ord(comp_char)]
            should_match = bool(comp_op(ref_char, query_char))
            if kmer_finder.kmers_present(comp_char) is not should_match:
                raise ValueError(
                    f"{char} should{' ' if should_match else ' not '}match {comp_char}"
                )


def test_kmer_finder_initialize_bigword():
    with pytest.raises(ValueError) as error:
        KmerFinder([(0, None, ["A" * 64])])
    error.match("A" * 64)
    error.match("64")


def test_kmer_finder_initialize_total_greater_than_max():
    kmer_finder = KmerFinder([(0, None, ["A" * 31, "B" * 31, "C" * 31])])
    assert kmer_finder.kmers_present("X" * 100 + "A" * 31)
    # assert kmer_finder.kmers_present("X" * 100 + "B" * 31)
    # assert kmer_finder.kmers_present("X" * 100 + "C" * 31)
    # assert not kmer_finder.kmers_present(string.ascii_letters)