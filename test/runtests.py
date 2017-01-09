"""
Unit tests of the shared library.
"""
from os import listdir
from os.path import isfile, join
import unittest
from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
import pairwise
# from Bio.SubsMat import MatrixInfo #  to be used later on

TEST_DATA_FOLDER = 'testdata'
DISPLAY_ALIGNMENT_LENGTH = 30
GAP_SYMBOL = '-'


def test_sequences_pairs():
    """
    Returns generator of all test data (as tuple: (filename, seq1, seq2)) that we have stored in testdata folder
    """
    seq_files = (f for f in listdir(TEST_DATA_FOLDER) if f.endswith('fasta') and isfile(join(TEST_DATA_FOLDER, f)))
    return (tuple([f] + [str(s.seq) for s in SeqIO.parse(join(TEST_DATA_FOLDER, f), 'fasta')]) for f in seq_files)


def shorten_string(string, max_length):
    """ Shorten string to given length adding ... if it is shorter"""
    if len(string) > max_length:
        return string[:max_length/2] + '...' + string[-max_length/2:]
    return string


class TestGlobalxx(unittest.TestCase):
    """
    Tests for globalxx method (no gap penalty, simply 1 point for match and 0 for mismatch)
    """

    def __init__(self, testname, filename, seq1, seq2):
        unittest.TestCase.__init__(self, testname)
        self.filename = filename
        self.seq1 = seq1
        self.seq2 = seq2

    def setUp(self):
        self.our_alignments = pairwise.globalxx(self.seq1, self.seq2)
        self.bio_alignments = pairwise2.align.globalxx(self.seq1, self.seq2)

    def print_alignment_error(self, err):
        """ Prints out info about assertion error about an alignment"""
        print 'Assertion error in ' + str(self.filename)
        print 'Error message: ' + str(err)
        print 'Seq1 = ' + shorten_string(self.seq1, 50)
        print 'Seq2 = ' + shorten_string(self.seq2, 50)

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


def print_alignments():
    """ Print all alignments. For debugging"""
    for filename, seq1, seq2 in test_sequences_pairs():
        our_alignments = pairwise.globalxx(seq1, seq2)
        bio_alignments = pairwise2.align.globalxx(seq1, seq2)

        print "File " + filename
        print "Our alignments:"
        print '\n'.join(format_alignment(*al) for al in our_alignments)
        print "Bio alignments:"
        print '\n'.join(format_alignment(*al) for al in bio_alignments)


def suite():
    """ Forms a test suite out of all test cases """
    t_suite = unittest.TestSuite()
    globalxx_tests = ['test_lengths', 'test_trailing_gaps', 'test_score_correctness', 'test_compare_bio_scores']
    for filename, seq1, seq2 in test_sequences_pairs():
        for testname in globalxx_tests:
            t_suite.addTest(TestGlobalxx(testname, filename, seq1, seq2))
    return t_suite


if __name__ == "__main__":
    # print_alignments()
    unittest.TextTestRunner().run(suite())
