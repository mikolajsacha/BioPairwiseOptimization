from timeit import default_timer as timer
from Bio import pairwise2
import pairwise
from unittests.alignment_test import test_sequences_pairs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_times(our_method, bio_method, seq1, seq2, *args, **kwargs):
    run_times = 5
    start = timer()
    for _ in xrange(run_times):
        our_method(seq1, seq2, *args, **kwargs)
    end = timer()
    our_time = ((end - start) * 1000) / run_times

    start = timer()
    for _ in xrange(run_times):
        bio_method(seq1, seq2, *args, **kwargs)
    end = timer()
    bio_time = ((end - start) * 1000) / run_times

    length = min(len(seq1), len(seq2))
    return length, our_time, bio_time


def run_test(description, our_method, bio_method, *args, **kwargs):
    results = [get_times(our_method, bio_method, s[1], s[2], *args, **kwargs) for s in test_sequences_pairs()]
    results.sort(key=lambda x: x[0])

    fig = plt.figure() 
    fig.canvas.set_window_title(description) 
    plt.title(description)
    plt.plot([r[0] for r in results], [r[1] for r in results], 'bo')
    plt.plot([r[0] for r in results], [r[2] for r in results], 'ro')
    plt.plot([r[0] for r in results], [r[1] for r in results], color='blue')
    plt.plot([r[0] for r in results], [r[2] for r in results], color='red')

    plt.ylabel('Execution time (ms)')
    plt.xlabel('Length of sequences')
    bio_plot_legend = mpatches.Patch(color='blue', label='Our performance')
    our_plot_legend = mpatches.Patch(color='red', label='Performance of Bio.pairwise2')
    plt.legend(handles=[our_plot_legend, bio_plot_legend])
    plt.show()


if __name__ == "__main__":
    print "Running globalxx test (scores only)..."
    run_test("globalxx (scores only)", pairwise.globalxx, pairwise2.align.globalxx, score_only=True)

    print "Running globalxx test (including alignments backtracking)..."
    run_test("globalxx (including alignments)", pairwise.globalxx, pairwise2.align.globalxx)
