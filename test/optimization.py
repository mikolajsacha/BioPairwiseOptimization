from timeit import default_timer as timer
from Bio import pairwise2
import pairwise
from unittests.alignment_test import test_sequences_pairs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
from itertools import cycle

PLOT_SAFE_FOLDER = "plots"

def fullname(o):
  return o.__module__ + "." + o.__name__

def get_times(methods, seq1, seq2, *args, **kwargs):
    run_times = 3
    times = []

    for method in methods:
        start = timer()
        for _ in xrange(run_times):
            method(seq1, seq2, *args, **kwargs)
        end = timer()
        time = ((end - start) * 1000) / run_times
        times.append(time)

    length = (len(seq1) + len(seq2)) / 2
    return length, times

def save_plot(directory, filename):
    if not os.path.exists(directory):
        os.makedirs(directory)
    plt.savefig(os.path.join(directory, filename))

def run_compare_test(description, compared_methods, *args, **kwargs):
    results = []
    for s in test_sequences_pairs():
        print("Testing {0}...".format(s[0]))
        result = get_times(compared_methods, s[1], s[2], *args, **kwargs)
        results.append(result)
    results.sort(key=lambda x: x[0])
    
    fig = plt.figure() 
    fig.canvas.set_window_title(description) 
    plt.title(description)
    colors = cycle(['r','g','b','yellow','magenta','cyan'])
    plot_legends = []

    for i, method in enumerate(compared_methods):
        color = next(colors)
        plt.plot([r[0] for r in results], [r[1][i] for r in results], 'bo', color=color)
        plt.plot([r[0] for r in results], [r[1][i] for r in results], color=color)
        plot_legends.append(mpatches.Patch(color=color, label=fullname(method)))

    plt.ylabel('Execution time (ms)')
    plt.xlabel('Length of sequences')
    plt.legend(handles=plot_legends)

    save_plot(PLOT_SAFE_FOLDER, description)
    print("Plot saved in \"{0}\"".format(PLOT_SAFE_FOLDER))

    plt.show()


if __name__ == "__main__":
    # print("Running globalxx test (scores only)...")
    # run_compare_test("globalxx (scores only)", [pairwise.globalxx, pairwise2.align.globalxx], score_only=True)

    # print("Running globalxx test (including alignments backtracking)...")
    # run_compare_test("globalxx (including alignments)", [pairwise.globalxx, pairwise2.align.globalxx])

    print("Test Bio.pairwise2 itself...")
    run_compare_test("Test Bio.pairwise2 itself...", [pairwise2.align.globalxx, pairwise2.align.localxx])
    run_compare_test("Test Bio.pairwise2 itself (score only)...", [pairwise2.align.globalxx, pairwise2.align.localxx], score_only=True)
