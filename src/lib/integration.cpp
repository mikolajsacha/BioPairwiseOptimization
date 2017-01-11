#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include <string>
#include "../include/model.h"
#include "../include/match_scorers.h"
#include "../include/global_alignment.h"

using namespace boost::python;
typedef std::vector<tuple> ATupleCollection;                      

/*! Converts an alignment to tuple - Bio.pairwise2 format */
tuple alignment_to_tuple(const Alignment& a)
{
    return make_tuple(a.sequence1, a.sequence2, a.score, a.begin, a.end);
}

/*! Converts alignments array to python list */
list align_array_to_py_list(Alignment* arr, unsigned length)
{
    list l;
    for (unsigned i = 0; i < length; i++) {
        l.append(alignment_to_tuple(arr[i]));
    }
    return l;
}

object globalxx(str seq1, str seq2, bool score_only=false) {
    std::string seq1_str = extract<std::string>(seq1);
    std::string seq2_str = extract<std::string>(seq2);
    ConstMatchScorer scorer(1,0);
    GlobalAlignment global_a(seq1_str, seq2_str);
    if (score_only) {
        global_a.populate_matrix_linear_gap_penalty_only_grid(&scorer, 0.0);
        return (object)global_a.get_score();
    }
    global_a.populate_matrix_linear_gap_penalty(&scorer, 0.0);
    global_a.backtrace_alignments();
    return align_array_to_py_list(global_a.alignments, global_a.alignments_count);
}

BOOST_PYTHON_MODULE(pairwise)
{
    class_<ATupleCollection >("ATupleCollection")
        .def(vector_indexing_suite<ATupleCollection >());
    def("globalxx", globalxx, (arg("seq1"), arg("seq2"), arg("score_only")=false));
}
