#include "match_scorers.h"
#include "global_alignment.h"
#include "model.h"
#include <iostream>
#include <vector>

int main() {
    MatchScorer * scorer = new ConstMatchScorer(1.0, -1.0);

    char const *seq1 = "ACTGTC";
    char const *seq2 = "ACGTGTC";

    std::vector<Alignment> alignments = global_alignment_linear_gap_penalty(seq1, 6, seq2, 7, scorer, -5);

    for (int i = 0; i < alignments.size(); i++) {
        std::cout<<alignments[i].sequence1<<std::endl;
        std::cout<<alignments[i].sequence2<<std::endl;
        std::cout<<std::endl;
    }
    std::cout<<"Alignment score: "<<alignments[0].score<<std::endl;

    delete scorer;
    return 0;
}
