#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H

#include "match_scorers.h"
#include "model.h"
#include <vector>

/*! Returns best global alignment score with a linear gap penalty and given score method */
/*!
  \param seq1 first sequence to align
  \param len1 length of the first sequence
  \param seq2 second sequence to align
  \param len2 length of the second sequence
  \param scorer a MatchScorer class
  \param penalty penalty for a gap
  \return vector of alignments with best score
*/
std::vector<Alignment> global_alignment_linear_gap_penalty(char const *seq1, unsigned len1, char const *seq2, unsigned len2, MatchScorer* scorer, double penalty);

#endif
