#!/usr/bin/env python3

import sys

import dnaio


def generate_kmers(sequence: str, k: int):
    for i in range(len(sequence) - k):
        yield sequence[i: i+k]


def main():
    adapter_kmers = set(generate_kmers(sys.argv[2], 8))
    with dnaio.open(sys.argv[1], mode="r", open_threads=0) as reader:
        for record in reader:
            sequence_kmers = set(generate_kmers(record.sequence, 8))
            if adapter_kmers.isdisjoint(sequence_kmers):
                pass


if __name__ == "__main__":
    main()