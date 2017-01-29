"""
Unit tests check if our optimized code gives exactly same results as original code.
It is a very important requirement to make our version of library is interchangeable with original.
"""
import unittest
import itertools
from timeit import default_timer as timer
from Bio import pairwise2
from lib import optimized_pairwise2 as opt_pairwise2
from testdata.test_sequences import get_test_sequences_pairs, all_equal


class AlignmentsEquivalenceTestCase(unittest.TestCase):
    """ Class with unit tests for testing if alignment methods give same results"""

    """ Names of all tested methods with tested arguments """
    test_cases = [
        ("globalxx", []),
        ("localxx", []),
        ("globalmx", [4, -2]),
        ("localmx", [4, -2]),
        ("globalms", [2, -1, -.5, -.1]),
        ("localms", [2, -1, -.5, -.1])
    ]

    """ Keyword arguments possible for all methods """
    keyword_arguments = {
        "penalize_extend_when_opening": [True, False],
        "penalize_end_gaps": [True, False],
        "one_alignment_only": [True, False],
    }

    def __init__(self, testname, seq1, seq2):
        unittest.TestCase.__init__(self, testname)
        self.seq1 = seq1
        self.seq2 = seq2

    def get_all_kwargs_combinations(self, score_only=False):
        params_values = [[(param, val) for val in values] for param, values in self.keyword_arguments.iteritems()]
        if score_only:
            for params in params_values:
                params.append(("score_only", True))
        return [dict(tuple_list) for tuple_list in itertools.product(*params_values)]

    def assert_alignments(self, results):
        """ Test that given lists of alignments are equal """

        # first assert that methods return same counts of alignments
        self.assertTrue(all_equal(list(map(len, results))))
        alignments_count = len(results[0])

        # test each returned alignment
        for i in xrange(alignments_count):
            # assert that alignment tuples have same lengths
            alignments = [r[i] for r in results]
            self.assertTrue(all_equal(list(map(len, alignments))))
            alignments_tuple_length = len(results[0][i])
            # assert that all alignments' tuples' elements are equal
            for j in xrange(alignments_tuple_length):
                alignments_fields = [al[j] for al in alignments]
                self.assertTrue(all_equal(alignments_fields))

    def run_tests(self, score_only=False):
        """ Runs tests for all parameters combinations """
        kwargs_combinations = self.get_all_kwargs_combinations(score_only=score_only)

        for method_name, args in AlignmentsEquivalenceTestCase.test_cases:
            methods = [pairwise2.align.__getattr__(method_name),
                       opt_pairwise2.align.__getattr__(method_name)]
            for kwargs in kwargs_combinations:
                # results - list of results from each method. Each result = list of alignments or score
                results = [method(self.seq1, self.seq2, *args, **kwargs) for method in methods]
                if score_only:
                    self.assertTrue(all_equal(results))
                else:
                    self.assert_alignments(results)


    def test_scores_equivalence(self):
        """ Test all test cases with all kwargs combinations for score with score_only=True option. """
        self.run_tests(score_only=True)

    def test_alignments_equivalence(self):
        """ Test all test cases with all kwargs combinations for equivalence of alignments. """
        self.run_tests()

    @staticmethod
    def get_test_suite(sequence_pairs):
        """ Forms a test suites for provided collection of samples """
        t_suite = unittest.TestSuite()
        for seq1, seq2 in sequence_pairs:
            t_suite.addTest(AlignmentsEquivalenceTestCase("test_scores_equivalence", seq1, seq2))
            t_suite.addTest(AlignmentsEquivalenceTestCase("test_alignments_equivalence", seq1, seq2))
        return t_suite

if __name__ == "__main__":
    print("Tested methods: ")
    print("Running unit tests (might take a while)...")
    print("\n".join("{0}, args: {1}".format(m, str(args)) for m, args in AlignmentsEquivalenceTestCase.test_cases))

    sequence_pairs = [(seq1, seq2) for _, seq1, seq2 in get_test_sequences_pairs(['unit'])]
    test_suite = AlignmentsEquivalenceTestCase.get_test_suite(sequence_pairs)

    unittest.TextTestRunner().run(test_suite)