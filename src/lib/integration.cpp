#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include <string>
#include "../include/model.h"
#include "../include/match_scorers.h"
#include "../include/global_alignment.h"

using namespace boost::python;
typedef std::vector<tuple> ATupleCollection;                      

/*! Converts std::vector to python list */
template<class T>
list std_vector_to_py_list(const std::vector<T>& v)
{
    object get_iter = iterator<std::vector<T> >();
    object iter = get_iter(v);
    list l(iter);
    return l;
}

/*! Converts an alignment to tuple - Bio.pairwise2 format */
tuple alignment_to_tuple(const Alignment& a)
{
    return make_tuple(a.sequence1, a.sequence2, a.score, a.begin, a.end);
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
    std::vector<Alignment> alignments = global_a.backtrace_alignments();
    std::vector<tuple> alignments_tuples;
    for (auto const &alignment : alignments) {
        alignments_tuples.push_back(alignment_to_tuple
                (alignment));
    }
    return std_vector_to_py_list(alignments_tuples);
}

BOOST_PYTHON_MODULE(pairwise)
{
    class_<ATupleCollection >("ATupleCollection")
        .def(vector_indexing_suite<ATupleCollection >());
    def("globalxx", globalxx, (arg("seq1"), arg("seq2"), arg("score_only")=false));
}
