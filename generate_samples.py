import random 
from numpy.random import normal
from Bio import pairwise2
import os.path

TEST_FOLDER = "testdata"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

def random_char_list(alphabet, length):
    return [random.choice(alphabet) for _ in xrange(length)]


def generate_sample(alphabet, length):
    seq1 = random_char_list(alphabet, length)
    seq2 = seq1[:]

    randomness = random.random()
    last_deleted = False

    i = 0
    while i < len(seq1) and i < len(seq2):
        if last_deleted and random.random() < 0.5:
            del seq2[i]
            continue
        last_deleted = False
        has_changed = random.random() < randomness
        if has_changed:
            change = random.randint(0,3)
            if change == 0: # change protein
                seq2[i] = random.choice(alphabet)
                i += 1
            if change == 1: # add protein
                seq2.insert(i, random.choice(alphabet))
                i += 2
            if change == 2: # remove protein
                del seq2[i]
                last_deleted = True
        else:
            i += 1

    added_front = int(normal(0, length**0.5))
    added_back = int(normal(0, length**0.5))

    if added_front > 0:
        seq2 = random_char_list(alphabet, added_front) + seq2
    elif added_front < 0:
        seq2 = seq2[-added_front:]

    if added_back > 0:
        seq2 += random_char_list(alphabet, added_front)
    elif added_back < 0:
        seq2 = seq2[:added_back]

    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2

def generate_test_samples(count, max_length):
    diff = max_length/count
    i = int(diff)
    while i <= max_length:
        yield generate_sample(PROTEIN_ALPHABET, i)
        i = int(i + diff)

def safe_samples(directory, samples):
    for seq1, seq2 in samples:
        seq_name = "generated_" + str(len(seq1))
        filename = seq_name + ".fasta"
        with open(os.path.join(directory, filename), 'w') as f:
            f.write(">{0}_seq1\n".format(seq_name))
            f.write(seq1 + '\n')
            f.write(">{0}_seq2\n".format(seq_name))
            f.write(seq2 + '\n')

if __name__ == "__main__":
    count = raw_input("Insert number of samples to generate: ")
    count = int(count)
    max_length = raw_input("Insert maximum length of a sample: ")
    max_length = int(max_length)
    samples = generate_test_samples(count, max_length)
    safe_samples(TEST_FOLDER, samples)


