"""Compute and compare bottom-s sketches of DNA sequences.
"""
import mmh3
from functools import partial
from sortedcontainers import SortedList

import dinopy as dp

mmhash = partial(mmh3.hash)


def bottom_sketch(seq, k, s, hf=mmhash):
    """Compute the bottom s sketch of the seqeunce using k-mers.

    Arguments:
        seq (bytes): Seqeunce to be sketched.
        k (int): k-mer size
        s (int): Sketch size, i.e. the number of hash values kept
                in the minimums list.
        hf (function): Hash function to be used. Default mmh3.hash
                with a fixed seed.

    Returns:
        SortedList: the s smallest hash values.
    """
    kmers = dp.qgrams(seq, k)
    # initialize the minimal hash value list with the first s values
    first_s_kmers = [kmers.__next__() for _ in range(s)]
    mins = SortedList([hf(kmer) for kmer in first_s_kmers])
    biggest_min = mins[-1]

    # traverse the remaining kmers
    for index, kmer in enumerate(kmers, s+1):
        hv = hf(kmer)

        if hv < biggest_min:
            # if the new has value is smaller than the biggest minimizer in
            # the list remove the biggest minimizer and add the new hv to
            # the sorted list
            mins.pop()  # remove the biggest minimizer (at the last position)
            mins.add(hv)
            biggest_min = mins[-1]

    return mins


def compare_sketches(sketch_a, sketch_b, sketch_size):
    """Estimate jaccard index by comparing number of shared sketch entries.
    """
    active_sketch_a = 0
    active_sketch_b = 0
    matches = 0
    for i in range(sketch_size):
        if sketch_a[active_sketch_a] == sketch_b[active_sketch_b]:
            matches += 1
            active_sketch_a += 1
            active_sketch_b += 1
        elif sketch_a[active_sketch_a] < sketch_b[active_sketch_b]:
            active_sketch_a += 1
        elif sketch_a[active_sketch_a] > sketch_b[active_sketch_b]:
            active_sketch_b += 1
    return matches / sketch_size
