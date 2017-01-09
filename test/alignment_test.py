"""
Basic class for unit testing alignments
"""
from os import listdir
from os.path import isfile, join
import unittest
from Bio import SeqIO
import pairwise

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


class AlignmentTest(unittest.TestCase):
    """
    Basic class for unit testing alignments
    """

    def __init__(self, testname, filename, seq1, seq2):
        unittest.TestCase.__init__(self, testname)
        self.filename = filename
        self.seq1 = seq1
        self.seq2 = seq2

    def print_alignment_error(self, err):
        """ Prints out info about assertion error about an alignment"""
        print 'Assertion error in ' + str(self.filename)
        print 'Error message: ' + str(err)
        print 'Seq1 = ' + shorten_string(self.seq1, 50)
        print 'Seq2 = ' + shorten_string(self.seq2, 50)

    @staticmethod
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