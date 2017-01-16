from timeit import default_timer as timer
from Bio.pairwise2 import align, affine_penalty, _find_start, _recover_alignments
import pairwise
from unittests.alignment_test import test_sequences_pairs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
from itertools import cycle
import optimization
import generate_samples
from Bio import BiopythonWarning 

# Try and load C implementations of functions. If I can't, 
# then throw a warning and use the pure Python implementations. 
try: 
    from cpairwise2 import rint, _make_score_matrix_fast 
except ImportError: 
    warnings.warn('Import of C module failed. Falling back to pure Python ' + 
                 'implementation. This may be slooow...', BiopythonWarning) 


class timed_align(align.alignment_function):
    preparation_time = 0
    make_score_matrix_time = 0
    find_starts_time = 0
    filter_starts_time = 0
    alignment_recovery_time = 0

    def __init__(self, name):
        align.alignment_function.__init__(self, name)

    @staticmethod
    def clear_timers():
        timed_align.preparation_time = 0
        timed_align.make_score_matrix_time = 0
        timed_align.find_starts_time = 0
        timed_align.filter_starts_time = 0
        timed_align.alignment_recovery_time = 0

    @staticmethod
    def align_timed(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
               penalize_extend_when_opening, penalize_end_gaps, 
               align_globally, gap_char, force_generic, score_only, 
               one_alignment_only): 
        """Return a list of alignments between two sequences or its score.
           Save times""" 

        # preparation -------------------------------------------------
        start = timer()

        if not sequenceA or not sequenceB: 
            return [] 
        try: 
            sequenceA + gap_char 
            sequenceB + gap_char 
        except TypeError: 
            raise TypeError('both sequences must be of the same type, either ' + 
                            'string/sequence object or list. Gap character must ' + 
                            'fit the sequence type (string or list)') 
     
        if not isinstance(sequenceA, list): 
            sequenceA = str(sequenceA) 
        if not isinstance(sequenceB, list): 
            sequenceB = str(sequenceB) 

        timed_align.preparation_time += (timer() - start) * 1000

        # making score matrix ----------------------------------------
        start = timer()
        if (not force_generic) and isinstance(gap_A_fn, affine_penalty) \
            and isinstance(gap_B_fn, affine_penalty): 
            open_A, extend_A = gap_A_fn.open, gap_A_fn.extend 
            open_B, extend_B = gap_B_fn.open, gap_B_fn.extend 
            x = _make_score_matrix_fast( 
                sequenceA, sequenceB, match_fn, open_A, extend_A, open_B, 
                extend_B, penalize_extend_when_opening, penalize_end_gaps, 
                align_globally, score_only) 
        else: 
            print("ATTENTION : Can't use C library for some reason.\nUsing slower pure Python version!")
            x = _make_score_matrix_generic( 
                sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
                penalize_end_gaps, align_globally, score_only) 
        score_matrix, trace_matrix = x 

        timed_align.make_score_matrix_time += (timer() - start) * 1000

        # finding starts ---------------------------------------------
        start = timer()
        # print("SCORE %s" % print_matrix(score_matrix)) 
        # print("TRACEBACK %s" % print_matrix(trace_matrix)) 
     
        # Look for the proper starting point. Get a list of all possible 
        # starting points. 
        starts = _find_start(score_matrix, align_globally) 

        timed_align.find_starts_time += (timer() - start) * 1000

        # filtering starts -------------------------------------------
        start = timer()

        # Find the highest score. 
        best_score = max([x[0] for x in starts]) 
     
        # If they only want the score, then return it. 
        if score_only: 
            return best_score 
     
        tolerance = 0  # XXX do anything with this? 
        # Now find all the positions within some tolerance of the best 
        # score. 
        starts = [(score, pos) for score, pos in starts 
                  if rint(abs(score - best_score)) <= rint(tolerance)] 
        timed_align.filter_starts_time += (timer() - start) * 1000

        # recovering alignments ---------------------------------------
        start = timer()
     
        # Recover the alignments and return them. 
        result = _recover_alignments(sequenceA, sequenceB, starts, score_matrix, 
                                    trace_matrix, align_globally, gap_char, 
                                    one_alignment_only, gap_A_fn, gap_B_fn) 
        timed_align.alignment_recovery_time += (timer() - start) * 1000
        return result

    @staticmethod
    def total_time():
        return timed_align.preparation_time + \
               timed_align.make_score_matrix_time + \
               timed_align.find_starts_time + \
               timed_align.filter_starts_time + \
               timed_align.alignment_recovery_time

    def __call__(self, *args, **keywds): 
                 keywds = self.decode(*args, **keywds) 
                 return timed_align.align_timed(**keywds) 



def draw_bar_chart(description):
    fig, ax = plt.subplots()
    values = [timed_align.preparation_time, timed_align.make_score_matrix_time,
              timed_align.find_starts_time, timed_align.filter_starts_time,
              timed_align.alignment_recovery_time]
    descriptions = ["Preparation", "Score matrix", "Find start", "Filter start", "Recover alignments"]
    width = 0.5
    rects = ax.bar(xrange(5), values, width, color='r', tick_label=descriptions)

    plt.title(description)
    plt.ylabel('Execution time (ms)')
    plt.xlabel('Bio.pairwise algorithm stage')

    optimization.save_plot(optimization.PLOT_SAFE_FOLDER, description)
    print("Plot saved in \"{0}\"".format(optimization.PLOT_SAFE_FOLDER))

    plt.show()

def perform_test(seq_len, seq_pairs, method_name, description, *args, **kwargs):
    print("-"*20)
    print(description)
    print("Using sequences of length ~ {0}".format(seq_len))
    timed_align.clear_timers()
    align = timed_align(method_name)
    for seq1, seq2 in seq_pairs:
        alignments = align(seq1, seq2, *args, **kwargs)

    seq_count = float(len(seq_pairs))
    timed_align.preparation_time /= seq_count
    timed_align.make_score_matrix_time /= seq_count
    timed_align.find_starts_time /= seq_count
    timed_align.filter_starts_time /= seq_count
    timed_align.alignment_recovery_time /= seq_count

    total_time = timed_align.total_time()
    print("-"*20)
    print("Total time: {:4.2f} ms".format(total_time))
    print("-"*20)
    print("Preparation: {:4.2f} ms".format(timed_align.preparation_time))
    print("Making score and trace matrices: {:4.2f} ms".format(timed_align.make_score_matrix_time))
    print("Finding possible starts: {:4.2f} ms".format(timed_align.find_starts_time))
    print("Filtering starts: {:4.2f} ms".format(timed_align.filter_starts_time))
    print("Recovering alignment: {:4.2f} ms".format(timed_align.alignment_recovery_time))
    draw_bar_chart(description)

if __name__ == "__main__":
    print("Testing how much time each part of algorithm takes in Bio.pairwise2...")
    seq_len = 3000
    samples_count = 3
    seq_pairs = [generate_samples.generate_sample(generate_samples.PROTEIN_ALPHABET, seq_len) for _ in xrange(samples_count)]

    perform_test(seq_len, seq_pairs, "localxx", "Analysis of Bio pairwise2 localxx execution time (score only)", score_only=True)
    perform_test(seq_len, seq_pairs, "globalxx", "Analysis of Bio pairwise2 globalxx execution time (score only)", score_only=True)
    perform_test(seq_len, seq_pairs, "localxx", "Analysis of Bio pairwise2 localxx execution time")
    perform_test(seq_len, seq_pairs, "globalxx", "Analysis of Bio pairwise2 globalxx execution time")
