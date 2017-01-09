from Bio import pairwise2
import pairwise
from alignment_test import test_sequences_pairs
from timeit import default_timer as timer
import matplotlib.pyplot as plt


def run_test(our_method, bio_method, seq1, seq2, *args, **kwargs):
    start = timer()
    our_method(seq1, seq2, *args, **kwargs)
    end = timer()
    our_time = (end - start) * 1000

    start = timer()
    bio_method(seq1, seq2, *args, **kwargs)
    end = timer()
    bio_time = (end - start) * 1000

    length = min(len(seq1), len(seq2))  # TODO is min or max better?
    return length, our_time, bio_time


def run_globalxx_test():
    results = [run_test(pairwise.globalxx, pairwise2.align.globalxx, s[1], s[2]) for s in test_sequences_pairs()]
    results.sort(key=lambda x: x[0])

    plt.plot([r[0] for r in results], [r[1] for r in results], color='blue')
    plt.plot([r[0] for r in results], [r[2] for r in results], color='red')

    plt.ylabel('Execution time (ms)')
    plt.xlabel('Length of sequences')
    plt.show()


if __name__ == "__main__":
    run_globalxx_test()
