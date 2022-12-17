#!/usr/bin/env python3

import argparse

import dnaio


from adapter_heuristic import create_kmers_and_offsets, kmers_present_in_sequence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adapter")
    parser.add_argument("fastq")
    parser.add_argument("with_adapter")
    parser.add_argument("no_adapter")
    args = parser.parse_args()
    kmers_and_offsets = create_kmers_and_offsets(args.adapter, 3, 0.1)
    print(kmers_and_offsets)
    with (
        dnaio.open(args.fastq, mode="r", open_threads=0) as reader,
        open(args.with_adapter, mode="wb") as with_adapter,
        open(args.no_adapter, mode="wb") as no_adapter,
    ):
        number_of_records = 0
        possible_adapters_found = 0
        for number_of_records, record in enumerate(reader, start=1):
            if kmers_present_in_sequence(kmers_and_offsets, record.sequence):
                with_adapter.write(record.fastq_bytes())
                possible_adapters_found += 1
            else:
                no_adapter.write(record.fastq_bytes())
    print(f"Percentage possible adapters: {possible_adapters_found / number_of_records}")


if __name__ == "__main__":
    main()