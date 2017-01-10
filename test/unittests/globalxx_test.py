"""
Tests for globalxx method (no gap penalty, simply 1 point for match and 0 for mismatch)
"""
import unittest
from Bio import pairwise2
import pairwise
from alignment_test import AlignmentTest, test_sequences_pairs, GAP_SYMBOL


class TestGlobalxx(AlignmentTest):
    """
    Tests for globalxx method (no gap penalty, simply 1 point for match and 0 for mismatch)
    """

    def __init__(self, testname, filename, seq1, seq2):
        AlignmentTest.__init__(self, testname, filename, seq1, seq2)

    def setUp(self):
        self.our_alignments = pairwise.globalxx(self.seq1, self.seq2)
        self.bio_alignments = pairwise2.align.globalxx(self.seq1, self.seq2)

    def test_lengths(self):
        """
        Tests if our alignments have correct lengths
        """
        try:
            for align in self.our_alignments:
                self.assertEqual(len(align[0]), len(align[1]))

        except AssertionError as err:
            print 'Testing globalxx lengths failed!'
            self.print_alignment_error(err)
            raise err

    def test_trailing_gaps(self):
        """
        Tests that alignments should not have trailing gaps ex: --AB-- and -AAB-- are incorrect
        """
        try:
            for align in self.our_alignments:
                self.assertFalse(align[0].startswith(GAP_SYMBOL) and align[0].endswith(GAP_SYMBOL)
                                 and align[1].startswith(GAP_SYMBOL) and align[1].endswith(GAP_SYMBOL))

        except AssertionError as err:
            print 'Testing globalxx trailing gaps failed!'
            self.print_alignment_error(err)
            raise err

    def test_score_correctness(self):
        """
        Check if scores for alignments are correct for given alignment sequences
        """
        try:
            for al1, al2, score, _, _ in self.our_alignments:
                calc_score = 0
                for i in xrange(len(al1)):
                    if al1[i] == al2[i]:
                        calc_score += 1

                self.assertEqual(score, calc_score)

        except AssertionError as err:
            print 'Testing scores correctness for globalxx failed!'
            self.print_alignment_error(err)
            raise err

    def test_compare_bio_scores(self):
        """
        Tests if scores of our alignments are same as scores calculated by methods from Bio.pairwise2
        """
        try:
            if self.bio_alignments and self.our_alignments:
                # assert that all our scores are same as scores from Bio.pairwise2
                for align in self.our_alignments:
                    self.assertEqual(align[2], self.bio_alignments[0][2])  # item no 2 is score

        except AssertionError as err:
            print 'Testing scores comparison with Bio.pairwise2 for globalxx failed!'
            self.print_alignment_error(err)
            raise err

    def test_compare_bio_lengths(self):
        """
        Our methods should not return less correct alignments than Bio.pairwise2 methods
        """
        try:
            # if bio_alignments are empty, assert that ours are also empty
            if not self.bio_alignments:
                self.assertFalse(self.our_alignments)
            else:
                self.assertGreaterEqual(len(self.our_alignments), len(self.bio_alignments))

        except AssertionError as err:
            print 'Testing alignments lengths comparison with Bio.pairwise2 for globalxx failed!'
            self.print_alignment_error(err)
            raise err

    @staticmethod
    def get_test_suite():
        """ Forms a test suite out of all test cases """
        t_suite = unittest.TestSuite()
        globalxx_tests = ['test_lengths', 'test_trailing_gaps', 'test_score_correctness', 'test_compare_bio_scores']
        for filename, seq1, seq2 in test_sequences_pairs(include_trivial=True):
            for testname in globalxx_tests:
                t_suite.addTest(TestGlobalxx(testname, filename, seq1, seq2))
        return t_suite

    @staticmethod
    def run_tests():
        test_suite = TestGlobalxx.get_test_suite()
        unittest.TextTestRunner().run(test_suite)

if __name__ == "__main__":
    # AlignmentTest.print_alignments()
    TestGlobalxx.run_tests()
