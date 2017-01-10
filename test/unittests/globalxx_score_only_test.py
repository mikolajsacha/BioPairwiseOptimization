"""
Tests for globalxx method (no gap penalty, simply 1 point for match and 0 for mismatch)
"""
import unittest
from Bio import pairwise2
import pairwise
from alignment_test import AlignmentTest, test_sequences_pairs, GAP_SYMBOL


class TestGlobalxxScoreOnly(AlignmentTest):
    """
    Tests for globalxx method with option score_only = True
    (no gap penalty, simply 1 point for match and 0 for mismatch)
    """

    def __init__(self, testname, filename, seq1, seq2):
        AlignmentTest.__init__(self, testname, filename, seq1, seq2)

    def setUp(self):
        self.our_score = pairwise.globalxx(self.seq1, self.seq2, score_only=True)
        self.bio_score = pairwise2.align.globalxx(self.seq1, self.seq2, score_only=True)

    def test_compare_bio_scores(self):
        """
        Tests if scores of our alignments are same as scores calculated by methods from Bio.pairwise2
        """
        try:
            # For empty dataset Bio.pairwise2 returns [] instead of some score, let's ignore this mismatch
            if self.seq1 and self.seq2:
                self.assertEqual(self.our_score, self.bio_score)
            else:
                self.assertEqual(self.our_score, 0.0)

        except AssertionError as err:
            print 'Testing scores comparison with Bio.pairwise2 for globalxx (score_only=True) failed!'
            self.print_alignment_error(err)
            raise err

    @staticmethod
    def get_test_suite():
        """ Forms a test suite out of all test cases """
        t_suite = unittest.TestSuite()
        globalxx_score_only_tests = ['test_compare_bio_scores']
        for filename, seq1, seq2 in test_sequences_pairs(include_trivial=True):
            for testname in globalxx_score_only_tests:
                t_suite.addTest(TestGlobalxxScoreOnly(testname, filename, seq1, seq2))
        return t_suite

    @staticmethod
    def run_tests():
        test_suite = TestGlobalxxScoreOnly.get_test_suite()
        unittest.TextTestRunner().run(test_suite)

if __name__ == "__main__":
    TestGlobalxxScoreOnly.run_tests()
