""" Analysis of original Bio.pairwise2 code performance """
from timeit import default_timer as timer
from Bio import pairwise2
import matplotlib.pyplot as plt
from Bio.pairwise2 import rint

import test_optimization as optimization
import generate_samples


class TimedAlign(pairwise2.align.alignment_function):
    """ Executes same alignment code as pairwise2.align, by measures time of its pairs for analysis """
    preparation_time = 0
    make_score_matrix_time = 0
    find_starts_time = 0
    filter_starts_time = 0
    alignment_recovery_time = 0

    def __init__(self, name):
        pairwise2.align.alignment_function.__init__(self, name)

    @staticmethod
    def align_timed(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
               penalize_extend_when_opening, penalize_end_gaps, 
               align_globally, gap_char, force_generic, score_only, 
               one_alignment_only): 
        """Return a list of alignments between two sequences or its score.
           This method is an exact copy of method from bio.pairwise2.align class.
           Only difference is that it measures times. """

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

        TimedAlign.preparation_time += (timer() - start) * 1000

        # making score matrix ----------------------------------------
        start = timer()
        if (not force_generic) and isinstance(gap_A_fn, pairwise2.affine_penalty) \
            and isinstance(gap_B_fn, pairwise2.affine_penalty): 
            open_A, extend_A = gap_A_fn.open, gap_A_fn.extend 
            open_B, extend_B = gap_B_fn.open, gap_B_fn.extend 
            x = pairwise2._make_score_matrix_fast( 
                sequenceA, sequenceB, match_fn, open_A, extend_A, open_B, 
                extend_B, penalize_extend_when_opening, penalize_end_gaps, 
                align_globally, score_only) 
        else: 
            print("ATTENTION : Can't use C library for some reason.\nUsing slower pure Python version!")
            x = pairwise2._make_score_matrix_generic( 
                sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
                penalize_end_gaps, align_globally, score_only) 
        score_matrix, trace_matrix = x 

        TimedAlign.make_score_matrix_time += (timer() - start) * 1000

        # finding starts ---------------------------------------------
        start = timer()
        # print("SCORE %s" % print_matrix(score_matrix)) 
        # print("TRACEBACK %s" % print_matrix(trace_matrix)) 
     
        # Look for the proper starting point. Get a list of all possible 
        # starting points. 
        starts = pairwise2._find_start(score_matrix, align_globally) 

        TimedAlign.find_starts_time += (timer() - start) * 1000

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
        TimedAlign.filter_starts_time += (timer() - start) * 1000

        # recovering alignments ---------------------------------------
        start = timer()
     
        # Recover the alignments and return them. 
        result = pairwise2._recover_alignments(sequenceA, sequenceB, starts, score_matrix, 
                                    trace_matrix, align_globally, gap_char, 
                                    one_alignment_only, gap_A_fn, gap_B_fn) 
        TimedAlign.alignment_recovery_time += (timer() - start) * 1000
        return result

    @staticmethod
    def clear_timers():
        """ Clears all timers for measuring execution times of algorithm steps """
        TimedAlign.preparation_time = 0
        TimedAlign.make_score_matrix_time = 0
        TimedAlign.find_starts_time = 0
        TimedAlign.filter_starts_time = 0
        TimedAlign.alignment_recovery_time = 0


    @staticmethod
    def total_time():
        """ Returns total measured execution time so far """
        return TimedAlign.preparation_time + \
               TimedAlign.make_score_matrix_time + \
               TimedAlign.find_starts_time + \
               TimedAlign.filter_starts_time + \
               TimedAlign.alignment_recovery_time

    def __call__(self, *args, **keywds):
        """ Overridden method from Bio.pairwise2 align class"""
        keywds = self.decode(*args, **keywds)
        return TimedAlign.align_timed(**keywds)



def draw_bar_chart(description):
    """ Plots bar chart of execution times measured so far by TimedAlign class """
    fig, ax = plt.subplots()
    values = [TimedAlign.preparation_time, TimedAlign.make_score_matrix_time,
              TimedAlign.find_starts_time, TimedAlign.filter_starts_time,
              TimedAlign.alignment_recovery_time]
    descriptions = ["Preparation", "Score matrix", "Find start", "Filter start", "Recover alignments"]
    width = 0.5
    ax.bar(xrange(5), values, width, color='r', tick_label=descriptions)

    plt.title(description)
    plt.ylabel('Execution time (ms)')
    plt.xlabel('Bio.pairwise algorithm stage')

    optimization.save_plot(optimization.PLOT_SAFE_FOLDER, description)
    print("Plot saved in \"{0}\"".format(optimization.PLOT_SAFE_FOLDER))

    plt.show()

def perform_test(seq_len, seq_pairs, method_name, description, *args, **kwargs):
    """ Performs full analysis of execution time for given sequence pairs"""
    print("-"*20)
    print(description)
    TimedAlign.clear_timers()
    align = TimedAlign(method_name)
    for seq1, seq2 in seq_pairs:
        alignments = align(seq1, seq2, *args, **kwargs)

    seq_count = float(len(seq_pairs))
    TimedAlign.preparation_time /= seq_count
    TimedAlign.make_score_matrix_time /= seq_count
    TimedAlign.find_starts_time /= seq_count
    TimedAlign.filter_starts_time /= seq_count
    TimedAlign.alignment_recovery_time /= seq_count

    total_time = TimedAlign.total_time()
    print("-"*20)
    print("Total time: {:4.2f} ms".format(total_time))
    print("-"*20)
    print("Preparation: {:4.2f} ms".format(TimedAlign.preparation_time))
    print("Making score and trace matrices: {:4.2f} ms".format(TimedAlign.make_score_matrix_time))
    print("Finding possible starts: {:4.2f} ms".format(TimedAlign.find_starts_time))
    print("Filtering starts: {:4.2f} ms".format(TimedAlign.filter_starts_time))
    print("Recovering alignment: {:4.2f} ms".format(TimedAlign.alignment_recovery_time))
    draw_bar_chart(description)
    print("")

if __name__ == "__main__":
    print("Testing how much time each part of algorithm takes in Bio.pairwise2...")
    seq_len = 3000
    print("Using sequence of length ~ {0}".format(seq_len))
    samples_count = 3
    seq_pairs = [generate_samples.generate_sample(generate_samples.PROTEIN_ALPHABET, seq_len) for _ in xrange(samples_count)]

    perform_test(seq_len, seq_pairs, "localxx", "Analysis of Bio pairwise2 localxx execution time (score only)", score_only=True)
    perform_test(seq_len, seq_pairs, "globalxx", "Analysis of Bio pairwise2 globalxx execution time (score only)", score_only=True)
    perform_test(seq_len, seq_pairs, "localxx", "Analysis of Bio pairwise2 localxx execution time")
    perform_test(seq_len, seq_pairs, "globalxx", "Analysis of Bio pairwise2 globalxx execution time")
