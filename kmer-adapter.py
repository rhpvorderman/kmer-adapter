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


def main():
    adapter_kmers = present_kmers(sys.argv[2], 0.1)
    print(adapter_kmers)
    with dnaio.open(sys.argv[1], mode="r", open_threads=0) as reader:
        number_of_records = 0
        possible_adapters_found = 0
        for number_of_records, record in enumerate(reader, start=1):
            sequence = record.sequence
            for kmer in adapter_kmers:
                if kmer in sequence:
                    possible_adapters_found += 1
                    break
    print(f"Percentage possible adapters: {possible_adapters_found / number_of_records}")



if __name__ == "__main__":
    main()