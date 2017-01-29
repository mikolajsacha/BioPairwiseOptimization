import os
from Bio import SeqIO

TEST_DATA_FOLDER = 'testdata'


def get_all_files(directory, extension):
    """ Returns generator of all files with provided extension from given directory """
    return (f for f in os.listdir(directory) if f.endswith(extension) and os.path.isfile(os.path.join(directory, f)))


def get_all_sequences(directory):
    """ Returns all test sequence pairs from given folder  (as tuples: (filename, seq1, seq2))"""
    files = get_all_files(directory, 'fasta')
    return (tuple([f] + [str(s.seq) for s in SeqIO.parse(os.path.join(directory, f), 'fasta')]) for f in files)


def all_equal(sequence):
    """ Returns true, if all elements of sequence are equal"""
    return all(x == sequence[0] for x in sequence)


def get_test_sequences_pairs(include_subfolders=[]):
    """ Returns list of all test data (as tuples: (filename, seq1, seq2))"""
    test_data = list(get_all_sequences(TEST_DATA_FOLDER))
    for folder in include_subfolders:
        folder_path = os.path.join(TEST_DATA_FOLDER, folder)
        test_data.extend(list(get_all_sequences(folder_path)))
    return test_data
