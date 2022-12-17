from typing import List, Tuple


def _chunk(sequence: str, chunks: int) -> List[str]:
    """Chunk a sequence in chunks parts"""
    chunk_size = len(sequence) // (chunks)
    remainder = len(sequence) % (chunks)
    offset = 0
    chunks = []
    for _ in range(remainder):
        chunks.append(sequence[offset:offset + chunk_size + 1])
        offset += (chunk_size + 1)
    while offset < len(sequence):
        chunks.append(sequence[offset:offset + chunk_size])
        offset += chunk_size
    return chunks


def create_kmers_and_offsets(adapter: str, min_overlap: int, error_rate: float
                             ) -> List[Tuple[str, int]]:
    adapter_length = len(adapter)
    max_errors = int(adapter_length * error_rate)
    error_lengths = []
    max_error = 1
    kmers_and_offsets = []
    for i in range(adapter_length + 1):
        if i * error_rate >= max_error:
            error_lengths.append(i)
            max_error += 1

    # Build up the array with chunks which should occur at the tail end
    # if the adapter overlaps with the end.
    min_overlap_kmer = adapter[:min_overlap]
    min_overlap_kmer_offset = -(error_lengths[0] - 1)
    kmers_and_offsets.append((min_overlap_kmer, min_overlap_kmer_offset))
    for i, error_length in enumerate(error_lengths):
        if (i + 1) < len(error_lengths):
            next_length = error_lengths[i + 1]
        else:
            next_length = adapter_length
        offset = -(next_length - 1)
        number_of_errors = i + 1
        chunks = _chunk(adapter[:error_length], number_of_errors + 1)
        chunk_offset = offset
        for chunk in chunks:
            kmers_and_offsets.append((chunk, chunk_offset))

    # Create kmers at least one of which should be in the read when there is a
    # an adapter
    chunks = _chunk(adapter, max_errors + 1)
    chunk_offset = 0
    for chunk in chunks:
        kmers_and_offsets.append((chunk, chunk_offset))
    return kmers_and_offsets


def kmers_present_in_sequence(kmers_and_offsets: List[Tuple[str, int]],
                              sequence: str):
    for kmer, offset in kmers_and_offsets:
        if sequence.find(kmer, offset) != -1:
            return True
    return False
