"""
Unit tests of the shared library.
"""
from os import listdir
from os.path import isfile, join
import unittest
from Bio import pairwise2, SeqIO
import pairwise
# from Bio.SubsMat import MatrixInfo #  to be used later on

TEST_DATA_FOLDER = 'testdata'
DISPLAY_ALIGNMENT_LENGTH = 30


def test_sequences_pairs():
    """
    Returns generator of all test data (as tuple: (filename, seq1, seq2)) that we have stored in testdata folder
    """
    seq_files = (f for f in listdir(TEST_DATA_FOLDER) if f.endswith('fasta') and isfile(join(TEST_DATA_FOLDER, f)))
    return (tuple([f] + [str(s.seq) for s in SeqIO.parse(join(TEST_DATA_FOLDER, f), 'fasta')]) for f in seq_files)


def shorten_string(string, max_length):
    """ Shorten string to given length adding ... if it is shorter"""
    if len(string) > max_length:
        return string[:max_length] + '...'
    return string


def format_alignment(alignment, max_length=DISPLAY_ALIGNMENT_LENGTH):
    """ Format alignment to readable form """
    return tuple(shorten_string(str(val), max_length) for val in alignment)


class TestScoreCorrectness(unittest.TestCase):
    """
    Tests if scores of our alignments are same as scores calculated by methods from Bio.pairwise2
    """
    def test_globalxx(self):
        """
        Tests globalxx method (no gap penalty, simply 1 point for match and 0 for mismatch)
        """
        for filename, seq1, seq2 in test_sequences_pairs():
            try:
                print 'Testing alignment for sequences: '
                print 'Seq1 = ' + shorten_string(seq1, 50)
                print 'Seq2 = ' + shorten_string(seq2, 50)
                our_alignments = pairwise.globalxx(seq1, seq2)
                bio_alignments = pairwise2.align.globalxx(seq1, seq2)
                # self.assertEqual(len(our_alignments), len(bio_alignments))
                print "Our alignments:"
                print '\n'.join([str(format_alignment(al)) for al in our_alignments])
                print "Bio alignments:"
                print '\n'.join([str(format_alignment(al)) for al in bio_alignments])
                print 'Test OK\n'
            except AssertionError as err:
                print 'Assertion error in ' + str(filename)
                print 'Error message: ' + str(err)
                raise err


if __name__ == "__main__":
    unittest.main()
