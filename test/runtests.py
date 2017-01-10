"""
Unit tests of the shared library.
"""
import unittest
from unittests.globalxx_test import TestGlobalxx
from unittests.globalxx_score_only_test import TestGlobalxxScoreOnly

if __name__ == "__main__":
    # AlignmentTest.print_alignments()
    TestGlobalxx.run_tests()
    TestGlobalxxScoreOnly.run_tests()
