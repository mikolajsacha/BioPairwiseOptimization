#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H

#include "match_scorers.h"
#include "model.h"
#include <vector>
#include <string>

/*! Contains methods for performing global alignment */
class GlobalAlignment {
    private:
        float** grid; /*! Contains scores */
        char** trace; /*! Contains directions of backtracing */
        const std::string seq1;
        const std::string seq2;

    public: 
        /*!
          \param seq1 first sequence to align
          \param seq2 second sequence to align
        */
        GlobalAlignment(const std::string& seq1, const std::string& seq2);
        ~GlobalAlignment();

        /*! Returns all alignments for already populated grid and trace matrices */
        /*!
          \return vector of alignments with best score
        */
        std::vector<Alignment> backtrace_alignments();

        /*! returns best alignment score (value in right-bottom corner of grid matrx) */
        float get_score();

        /* runs Needleman - Wunsch algorithm, but populates only grid matrix.
         * Useful for score_only=True */
        /*!
          \param scorer a MatchScorer object
          \param penalty a linear gap penalty
        */
        void populate_matrix_linear_gap_penalty_only_grid(MatchScorer* scorer, float penalty);

        /*! runs dynamic Needleman - Wunsch algorithm populating grid and trace matrices */
        /*!
          \param scorer a MatchScorer object
          \param penalty a linear gap penalty
        */
        void populate_matrix_linear_gap_penalty(MatchScorer* scorer, float penalty);
};
#endif
