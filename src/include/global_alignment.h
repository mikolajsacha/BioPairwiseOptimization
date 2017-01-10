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

        /*! Returns all alignments for already populated grid and dir matrices */
        /*!
          \return vector of alignments with best score
        */
        std::vector<Alignment> backtrace_alignments();

        /*! returns best alignment score (value in right-bottom corner of grid matrx) */
        double get_score();

        /*! runs dynamic Needleman - Wunsch algorithm populating grid and dir matrices */
        /*!
          \param scorer a MatchScorer object
          \param penalty a linear gap penalty
        */
        void populate_matrix_linear_gap_penalty(MatchScorer* scorer, double penalty);
};


#endif
