# cython: language_level=3

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy, strstr, memset
from cpython.unicode cimport PyUnicode_CheckExact, PyUnicode_GET_LENGTH

cdef extern from "Python.h":
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)

ctypedef struct KmerEntry:
    size_t kmer_offset
    size_t kmer_length
    ssize_t search_offset


cdef class KmerFinder:
    cdef:
        char *kmers
        KmerEntry *kmer_entries
        size_t number_of_kmers

    def __cinit__(self, kmers_and_offsets):
        self.kmers = NULL
        self.kmer_entries = NULL
        self.number_of_kmers = 0
        kmers = [kmer for kmer, _ in kmers_and_offsets]
        kmer_total_length = sum(len(kmer) for kmer in kmers)
        number_of_entries = len(kmers_and_offsets)
        self.kmer_entries = <KmerEntry *>PyMem_Malloc(number_of_entries * sizeof(KmerEntry))
        # for the kmers the NULL bytes also need space.
        self.kmers = <char *>PyMem_Malloc(kmer_total_length + number_of_entries)
        self.number_of_kmers = number_of_entries
        cdef size_t kmer_offset = 0
        cdef char *kmer_ptr
        cdef Py_ssize_t kmer_length
        for i, (kmer, offset) in enumerate(kmers_and_offsets):
            if not PyUnicode_CheckExact(kmer):
                raise TypeError(f"Kmer should be a string not {type(kmer)}")
            if not PyUnicode_IS_COMPACT_ASCII(kmer):
                raise ValueError("Only ASCII strings are supported")
            self.kmer_entries[i].kmer_offset = kmer_offset
            self.kmer_entries[i].search_offset  = offset
            kmer_length = PyUnicode_GET_LENGTH(kmer)
            self.kmer_entries[i].kmer_length = kmer_length;
            kmer_ptr = <char *>PyUnicode_DATA(kmer)
            memcpy(self.kmers + kmer_offset, kmer_ptr, kmer_length)
            kmer_offset += kmer_length
            self.kmers[kmer_offset] = 0
            kmer_offset += 1
        print(self.kmers[0:kmer_total_length + number_of_entries])

    def kmers_present(self, str sequence):
        cdef:
            KmerEntry entry
            size_t i
            size_t kmer_offset
            size_t kmer_length
            ssize_t search_offset
            char *kmer_ptr
            char *search_ptr
            char *search_result
            size_t search_length
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError("Only ASCII strings are supported")
        cdef char *seq = <char *>PyUnicode_DATA(sequence)
        cdef Py_ssize_t seq_length = PyUnicode_GET_LENGTH(sequence)
        for i in range(self.number_of_kmers):
            entry = self.kmer_entries[i]
            search_offset = entry.search_offset
            if search_offset < 0:
                search_offset = seq_length + search_offset
                if search_offset < 0:
                    search_offset = 0
            if search_offset > seq_length:
                continue
            kmer_length = entry.kmer_length
            kmer_offset = entry.kmer_offset
            kmer_ptr = self.kmers + kmer_offset
            search_ptr = seq + search_offset
            search_length = seq_length - (search_ptr - seq)
            search_result = bitap_bitwise_search(search_ptr, search_length,
                                                 kmer_ptr, kmer_length)
            if search_result:
                return True
        return False

    def __dealloc__(self):
        PyMem_Free(self.kmers)
        PyMem_Free(self.kmer_entries)


cdef populate_needle_bitmap(size_t needle_bitmap[128], char *needle, size_t needle_length):
    cdef size_t i
    memset(needle_bitmap, 0xff, sizeof(size_t) * 128)
    for i in range(needle_length):
        needle_bitmap[needle[i]] &= ~(1UL << i)


cdef char *bitap_bitwise_search(char *haystack, size_t haystack_length,
                                char *needle, size_t needle_length):
    cdef:
        size_t R
        size_t pattern_mask[128]
        size_t i

    if needle_length == 0:
        return haystack
    if needle_length > (sizeof(size_t) * 8 -1 ):
        return "The pattern is too long!"

    # Initialize the bit array R
    R = ~1

    # Initialize the pattern bitmasks
    populate_needle_bitmap(pattern_mask, needle, needle_length)

    for i in range(haystack_length):
        # Update the bit array
        R |= pattern_mask[haystack[i]]
        R <<= 1

        if (0 == (R & (1UL << needle_length))):
            return (haystack + i - needle_length) + 1;

    return NULL