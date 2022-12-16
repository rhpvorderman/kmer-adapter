#!/usr/bin/env python3

import sys

import dnaio


def generate_kmers(sequence: str, k: int):
    for i in range(len(sequence) - k):
        yield sequence[i: i+k]


def present_kmers(sequence, error_rate):
    length = len(sequence)
    max_errors = int(length * error_rate)
    k = length // (max_errors + 1)
    return [sequence[i:i+k] for i in  range(0, length, k)][:max_errors + 1]


def illumina_truseq_candidate(sequence):
    """
    Hard-coded heuristic for the illumina TruSeq adapter.
    :param sequence: The sequence to check
    :return: Whether the TruSeq adapter might align.
    """
    # Values are hard-coded for speed and because this is a test
    # No errors: 1-9
    # Max 1 errors: 10-19
    # Max 2 errors: 20-29
    # Max 3 errors: 30-33
    # Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    # With min overlap 3 and no errors in 9 bases, if the adapter aligns at
    # the very end the sequence "AGA" should be present in the last 9 bases.
    if sequence.find("AGA", -9):
        return True
    # if more than 10 bases overlap, there can be one error. this means either
    # AGATC or GGAAG is present in the last 19 bases
    if sequence.find("AGATC", -19) or sequence.find("GGAAG", -19):
        return True
    # If more than 20 bases overlap either "AGATCGG", "AAGAGCA" or "CACGTC"
    # should be present, as there may be two errors.
    if (
            # First two sequences disabled because of overlap with next clause
            #sequence.find("AGATCGG", -29) or
            #sequence.find("AAGAGCA", -29) or
            sequence.find("CACGTC", -29)
    ):
        return True
    # If there are more than 30 bases overlap there may be 3 errors. So there
    # can be 4 pieces: AGATCGG, AAGAGCA, CACGTCTG, AACTCCAG
    if sequence.find("AGATCGG", -33) or sequence.find("AAGAGCA", -33) or sequence.find("CACGTCTG", -33) or  sequence.find("AACTCCAG", -33):
        return True
    # If the entire adapter is present one of the following 8-mers should be present:
    # AGATCGGA, AGAGCACA, CGTCTGAA, CTCCAGTC
    if sequence.find("AGATCGGA") or sequence.find("AGAGCACA") or sequence.find("CGTCTGAA") or sequence.find("CTCCAGTC"):
        return True
    return False


def main():
    adapter_kmers = present_kmers(sys.argv[2], 0.1)
    print(adapter_kmers)
    with dnaio.open(sys.argv[1], mode="r", open_threads=0) as reader:
        number_of_records = 0
        possible_adapters_found = 0
        for number_of_records, record in enumerate(reader, start=1):
            if illumina_truseq_candidate(record.sequence):
                possible_adapters_found += 1
    print(f"Percentage possible adapters: {possible_adapters_found / number_of_records}")


if __name__ == "__main__":
    main()