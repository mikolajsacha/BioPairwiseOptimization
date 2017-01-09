"""
Unit tests of the shared library.
"""
import unittest
from globalxx_test import TestGlobalxx

if __name__ == "__main__":
    # AlignmentTest.print_alignments()
    test_suite = TestGlobalxx.get_test_suite()
    unittest.TextTestRunner().run(test_suite)
