#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
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

list globalxx(str seq1_obj, str seq2_obj) {
    char const* seq1 = extract<const char*>(seq1_obj);
    char const* seq2 = extract<const char*>(seq1_obj);
    ConstMatchScorer scorer(1,0);
    std::vector<Alignment> alignments = global_alignment_linear_gap_penalty(seq1, len(seq1_obj), seq2, len(seq2_obj), &scorer, 0.0);
    std::vector<tuple> alignments_tuples;
    for (auto const &alignment : alignments) {
        alignments_tuples.push_back(alignment_to_tuple
                (alignment));
    }
    return std_vector_to_py_list(alignments_tuples);
}

BOOST_PYTHON_MODULE(pairwise)
{
    // class_<Alignment>("Alignment")
        // .def_readwrite("sequence1", &Alignment::sequence1)
        // .def_readwrite("sequence2", &Alignment::sequence2)
        // .def_readwrite("score", &Alignment::score)
        // .def_readwrite("begin", &Alignment::begin)
        //
    class_<ATupleCollection >("ATupleCollection")
        .def(vector_indexing_suite<ATupleCollection >());
    def("globalxx", globalxx);
}
