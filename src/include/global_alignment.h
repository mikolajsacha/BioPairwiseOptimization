#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H

#include "match_scorers.h"
#include "model.h"
#include <vector>
#include <string>

/*! Contains methods for performing global alignment */
class GlobalAlignment {
    private:
        std::vector<std::vector<double> > grid; /*! Contains scores */
        std::vector<std::vector<AlignmentDirection> > dir; /*! Contains directions of backtracing */
        std::string seq1;
        std::string seq2;

    public: 
        /*!
          \param seq1 first sequence to align
          \param seq2 second sequence to align
        */
        GlobalAlignment(std::string seq1, std::string seq2);

        /*! Returns best global alignment score with a linear gap penalty and given score method */
        /*!
          \param scorer a MatchScorer class
          \param penalty penalty for a gap
          \return vector of alignments with best score
        */
        std::vector<Alignment> align_linear_gap_penalty(MatchScorer* scorer, double penalty);
};


#endif
