# BioPairwiseOptimization

This is an educational project which aims to rewrite functionalities from pairwise2 module from BioPython library. Pairwise2 module provides methods for finding global and local allignments in string sequences which are commonly used in bioinformatics. Idea for optimizing this module appeared after discovering that it is often considerably slow for big inputs.

My planned method as of now is to implement internal algorithms in C++, using multithreading if possible, to create fastest code possible. Then C++ code will be integrated into Python. In Python, same interface as in original Bio.pairwise2 module should be possible to use, so switching to using my library is painless.
