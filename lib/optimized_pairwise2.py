import itertools
import Bio.pairwise2 as pairwise2

def _find_start(score_matrix, best_score): 
    """Return a list of starting points (score, (row, col)). 
    Indicating every possible place to start the tracebacks. 
    """ 
    nrows, ncols = len(score_matrix), len(score_matrix[0]) 

    return [(best_score, (row, col)) for row, col in itertools.product(xrange(nrows), xrange(ncols))
            if score_matrix[row][col] == best_score]

def _align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
           penalize_extend_when_opening, penalize_end_gaps, 
           align_globally, gap_char, force_generic, score_only, 
           one_alignment_only):
    """Return a list of alignments between two sequences or its score.
       Use optimized methods where possible""" 

    if not sequenceA or not sequenceB: 
        return [] 
    try: 
        sequenceA + gap_char 
        sequenceB + gap_char 
    except TypeError: 
        raise TypeError('both sequences must be of the same type, either ' + 
                        'string/sequence object or list. Gap character must ' + 
                        'fit the sequence type (string or list)') 
 
    if not isinstance(sequenceA, list): 
        sequenceA = str(sequenceA) 
    if not isinstance(sequenceB, list): 
        sequenceB = str(sequenceB) 

    if (not force_generic) and isinstance(gap_A_fn, pairwise2.affine_penalty) \
        and isinstance(gap_B_fn, pairwise2.affine_penalty): 
        open_A, extend_A = gap_A_fn.open, gap_A_fn.extend 
        open_B, extend_B = gap_B_fn.open, gap_B_fn.extend 
        x = pairwise2._make_score_matrix_fast( 
            sequenceA, sequenceB, match_fn, open_A, extend_A, open_B, 
            extend_B, penalize_extend_when_opening, penalize_end_gaps, 
            align_globally, score_only) 
    else: 
        x = pairwise2._make_score_matrix_generic( 
            sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, 
            penalize_end_gaps, align_globally, score_only) 
    score_matrix, trace_matrix = x 

    # Find the highest score. 
    nrows, ncols = len(score_matrix), len(score_matrix[0]) 
 
    if align_globally:
        starts = [(score_matrix[nrows-1][ncols-1], (nrows-1,ncols-1))]
        if score_only: 
            return starts[0][0]
    else:
        best_score = max(score_matrix[row][col] for row, col in itertools.product(xrange(nrows), xrange(ncols)))
        if score_only: 
            return best_score 

        starts = _find_start(score_matrix, best_score) 

    return pairwise2._recover_alignments(sequenceA, sequenceB, starts, score_matrix, 
                                         trace_matrix, align_globally, gap_char, 
                                         one_alignment_only, gap_A_fn, gap_B_fn) 

class align(object):
    """This class provides same functionalities as Bio.pairwise2 align class, but with some methods overloaded as optimized.""" 

    class alignment_function(pairwise2.align.alignment_function):
        def __init__(self, name):
            pairwise2.align.alignment_function.__init__(self, name)

        def __call__(self, *args, **keywds): 
                     keywds = self.decode(*args, **keywds) 
                     return _align(**keywds) 

    def __getattr__(self, attr): 
        return self.alignment_function(attr) 

align = align()
