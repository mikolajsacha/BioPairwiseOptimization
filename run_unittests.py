"""
Unit tests of the shared library.
"""
from abc import ABCMeta
import unittest
from Bio import pairwise2
from lib import optimized_pairwise2 as opt_pairwise2
from testdata.test_sequences import get_test_sequences_pairs, all_equal


class BaseAlignmentTestCase(unittest.TestCase):
    """ Abstract class for unit tests of alignment methods"""
    __metaclass__ = ABCMeta

    def __init__(self, testname, seq1, seq2):
        unittest.TestCase.__init__(self, testname)
        self.seq1 = seq1
        self.seq2 = seq2

    def test_methods_alignments(self, methods, *args, **kwargs):
        """ Test that given methods return same alignments """
        results = [method(*args, **kwargs) for method in methods]
        
        # first assert that methods return same counts of alignments
        self.assertTrue(all_equal(list(map(len, results))))
        alignments_count = len(results[0])

        for i in xrange(alignments_count):
            # assert that alignment tuples have same lengths
            alignments = [r[i] for r in results]
            self.assertTrue(all_equal(list(map(len, alignments))))
            alignments_tuple_length = len(results[0][i])
            for j in xrange(alignments_tuple_length):
                alignments_fields = [al[j] for al in alignments]
                self.assertTrue(all_equal(alignments_fields))

    @staticmethod
    def get_test_suite_for_tests(test_case_class, tests):
        """ Forms a test suite out of provided tests """
        t_suite = unittest.TestSuite()
        for _, seq1, seq2 in get_test_sequences_pairs(['unit']):
            for testname in tests:
                t_suite.addTest(test_case_class(testname, seq1, seq2))
        return t_suite

class GlobalxxTest(BaseAlignmentTestCase):
    """ Test globalxx method """

    def __init__(self, testname, seq1, seq2):
        BaseAlignmentTestCase.__init__(self, testname, seq1, seq2)

    def test_globalxx(self):
        """ Test if our and original globalxx give same results"""
        methods = [pairwise2.align.globalxx, opt_pairwise2.align.globalxx]
        self.test_methods_alignments(methods, self.seq1, self.seq2)

    @staticmethod
    def get_test_suite():
        return BaseAlignmentTestCase.get_test_suite_for_tests(
                GlobalxxTest, ['test_globalxx']
                )

class LocalxxTest(BaseAlignmentTestCase):
    """ Test localxx method """

    def __init__(self, testname, seq1, seq2):
        BaseAlignmentTestCase.__init__(self, testname, seq1, seq2)

    def test_localxx(self):
        """ Test if our and original localxx give same results"""
        methods = [pairwise2.align.localxx, opt_pairwise2.align.localxx]
        self.test_methods_alignments(methods, self.seq1, self.seq2)

    @staticmethod
    def get_test_suite():
        return BaseAlignmentTestCase.get_test_suite_for_tests(
                LocalxxTest, ['test_localxx']
                )

if __name__ == "__main__":
    testCases = [GlobalxxTest, LocalxxTest]

    for testCaseClass in testCases:
        test_suite = testCaseClass.get_test_suite()
        unittest.TextTestRunner().run(test_suite)
