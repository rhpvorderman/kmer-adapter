# from libc.stdint cimport uint16_t
#
#
#
# ctypedef struct _CheckMer:
#     char *sequence
#     Py_ssize_t offset


class AdapterHeuristic:
    # cdef:
    #     object sequences_obj
    #     char *sequences
    #     _CheckMer *checkmers

    @staticmethod
    def _chunk(sequence, chunks):
        """Chunk a sequence in chunks parts"""
        chunk_size = len(sequence) // (chunks)
        remainder = len(sequence) % (chunks)
        offset = 0
        chunks = []
        for _ in range(remainder):
            chunks.append(sequence[offset:offset+chunk_size+1])
            offset += (chunk_size + 1)
        while offset < len(sequence):
            chunks.append(sequence[offset:offset+chunk_size])
            offset += chunk_size
        return chunks

    def __init__(self, adapter, min_overlap, error_rate):
        adapter_length = len(adapter)
        max_errors = int(adapter_length * error_rate)
        error_lengths = []
        max_error = 1
        offsets_and_kmers = []
        for i in range(adapter_length + 1):
            if i * error_rate >= max_error:
                error_lengths.append(i)
                max_error += 1

        # Build up the array with chunks which should occur at the tail end
        # if the adapter overlaps with the end.
        min_overlap_kmer = adapter[:min_overlap]
        min_overlap_kmer_offset = -(error_lengths[0] -1)
        offsets_and_kmers.append((min_overlap_kmer_offset, min_overlap_kmer))
        for i, error_length in enumerate(error_lengths):
            if (i + 1) < len(error_lengths):
                next_length = error_lengths[i + 1]
            else:
                next_length = adapter_length
            offset = -(next_length - 1)
            chunks = self._chunk(adapter[:error_length], i + 2)
            chunk_offset = offset
            for chunk in chunks:
                offsets_and_kmers.append((chunk_offset, chunk))

        # Create kmers at least one of which should be in the read when there is a
        # an adapter
        chunks = self._chunk(adapter, max_errors + 1)
        chunk_offset = 0
        for chunk in chunks:
            offsets_and_kmers.append((chunk_offset, chunk))
        self.offsets_and_kmers = offsets_and_kmers

    def adapter_possibly_present(self, sequence):
        for offset, kmer in self.offsets_and_kmers:
            if sequence.find(kmer, offset) != -1:
                return True
        return False
