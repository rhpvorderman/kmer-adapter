# kmer-adapter

This repository was created to prove a concept for quicker adapter matching
in cutadapt. The current iteration of cutadapt:

- Allows a 0.1 error rate when matching the adapter.
- Takes into account indels

To do so effectively, it aligns the adapter to the sequence using a custom
semi-global alignment algorithm which runs in O(mn) time where m is the size
of the adapter and n is the size of the sequence.

When reading something about k-mers, it occurred to me that these could be used
to prove an adapter was **not** present in a much faster fashion. Given an
adapter with length 33 (TruSeq Illumina adapter), 3 errors are allowed in
cutadapt. When this adapter is cut in 4 (3+1) non-overlapping pieces, it
follows that at least one of the 4 pieces must occur without error. Since
modern string searching algorithms can search in approximately O(n) time this
reduces the calculation time to O((x+1)n) where x is the number of errors
allowed. Or O((e*m + 1)n) where e is the error rate. This is faster than
O(mn).

## Steps

- [x] make a proof of concept where the heuristic is proven using Python's 
  `.find()` method on strings.
- [x] cythonize the algorithm. All the kmers are saved in their C-form and
  passed to libc's `strstr` in order to remove Python overhead. This made
  the heuristic 3 times faster.
- String search algorithms that operate in O(n) time usually initialize some
  sort of search structure to allow this. Since the same sequences are sought 
  over and over again this search structure can be build up only once. Search
  an easily implementable search algorithm. Found the shift-OR algorithm 
  as a very fast candidate in this paper: https://arxiv.org/pdf/1012.2547v1.pdf
- [x] implement shift OR algorithm with the implementation on wikipedia:
  https://en.wikipedia.org/wiki/Bitap_algorithm.
- [x] Implement the initialization of the search structures (arrays of bitmasks)
  only once rather than for every string.
- [x] Shift-OR algorithm allows case-independent matching and IUPAC matching at
  no cost. Implement this.
- Try various optimizations:
  - [x] branchless programming. It is possible to return a boolean rather than 
    a char pointer since only the presence or non-presence is important. It 
    is therefore possible to bitwise AND everything together and only check at
    the end if it is zero. This eliminates an if-comparison in the loop. 
    This is **slower**. The slowest part in the loop is the bitmask lookup
    due to the lookup table being ASCII_WIDTH * sizeof(size_t) big (128 * 8) == 1KB.
    This is larger than a cacheline (64 bytes) but smaller than the l1 cache. 
    Expected lookup time is therefore 3-4 clockticks. This is slower than an
    if-comparison. Therefore for every item in the loop and exiting early is
    faster.
  - [x] Different integer sizes. In theory, we can fit 32 16-bit integers on 
    a 64 byte cacheline. In theory we can therefore do the shift-or search
    with the same cacheline if all characters compared are ACGT only. In 
    practice `size_t` integers (64-bit on 64-bit architecture) is the fastest. 
    Probably because the cache-line advantage is not there and the bitwise
    operations are faster using the native register size.
  