#include "global_alignment.h"
#include "match_scorers.h"
#include "model.h"
#include <vector>
#include <algorithm>

std::vector<Alignment> global_alignment_linear_gap_penalty(char const * seq1, int len1, char const * seq2, int len2, MatchScorer* scorer, double penalty) {
    double** grid = new double*[len1+1];
    AlignmentDirection** dir = new AlignmentDirection*[len1];

    for (int i = 0; i < len1; i++) {
        grid[i] = new double[len2+1];
        dir[i] = new AlignmentDirection[len2];
    }
    grid[len1] = new double[len2+1];

    for (int i = 0; i <= len1; i++) {
        grid[i][0] = i * penalty;
    }
    for (int i = 0; i <= len2; i++) {
        grid[0][i] = i * penalty;
    }

    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            double matchScore = scorer->getScore(seq1[i-1], seq2[j-1]);
            double diagScore = grid[i-1][j-1] + matchScore;
            double upScore = grid[i-1][j] + penalty;
            double leftScore = grid[i][j-1] + penalty;

            if (diagScore > upScore) {
                if (diagScore > leftScore) { //upScore < leftScore < diagScore
                    grid[i][j] = diagScore;
                    dir[i-1][j-1] = diagAlign;
                }
                else if (diagScore == leftScore) { // upScore < diagScore == leftScore
                    grid[i][j] = diagScore;
                    dir[i-1][j-1] = diagLeftAlign;
                }
                else { // upScore < diagScore < leftScore
                    grid[i][j] = leftScore;
                    dir[i-1][j-1] = leftAlign;
                }
            }
            else if (diagScore == upScore) {
                if (diagScore > leftScore) { //leftScore < upScore == diagScore
                    grid[i][j] = diagScore;
                    dir[i-1][j-1] = diagUpAlign;
                }
                else if (diagScore == leftScore) { // upScore == diagScore == leftScore
                    grid[i][j] = diagScore;
                    dir[i-1][j-1] = allAlign;
                }
                else { // upScore == diagScore < leftScore
                    grid[i][j] = leftScore;
                    dir[i-1][j-1] = leftAlign;
                }
            }
            else { // diagScore < upScore
                if (upScore > leftScore) { // diagScore < upScore, leftScore < upScore
                    grid[i][j] = upScore;
                    dir[i-1][j-1] = upAlign;
                }
                else if (upScore == leftScore) { // diagScore < upScore == leftScore
                    grid[i][j] = upScore;
                    dir[i-1][j-1] = upLeftAlign;
                }
                else { // diagScore < upScore < leftScore
                    grid[i][j] = leftScore;
                    dir[i-1][j-1] = leftAlign;
                }
            }
        }
    }

    /* Tracking back best alignments */
    std::vector<AlignmentBacktrace> backtraces;

    AlignmentBacktrace initBacktrace;
    initBacktrace.alignment1 = "";
    initBacktrace.alignment2 = "";
    initBacktrace.i1 = len1-1;
    initBacktrace.i2 = len2-1;
    initBacktrace.isFinished = false;
    backtraces.push_back(initBacktrace);

    int finished_backtraces = 0;
    while (finished_backtraces < backtraces.size()) {
        int backtracesSize = backtraces.size();
        for (int i = 0; i < backtracesSize; i++) {
            if (backtraces[i].isFinished) {
                continue;
            }
            if (backtraces[i].i1 < 0 || backtraces[i].i2 < 0) {
                backtraces[i].isFinished = true;
                finished_backtraces++;
                continue;
            }
            int i1 = backtraces[i].i1;
            int i2 = backtraces[i].i2;
            AlignmentDirection alignDir = dir[i1][i2];

            if ((alignDir & leftAlign) == leftAlign) {
                backtraces[i].alignment1.push_back('-');
                backtraces[i].alignment2.push_back(seq2[i2]);
                if ((alignDir & upAlign) == upAlign) {
                    AlignmentBacktrace newTrace(backtraces[i]);
                    newTrace.alignment1.push_back(seq1[i1]);
                    newTrace.alignment2.push_back('-');
                    newTrace.i1--;
                    backtraces.push_back(newTrace);
                }
                if ((alignDir & diagAlign) == diagAlign) {
                    AlignmentBacktrace newTrace(backtraces[i]);
                    newTrace.alignment1.push_back(seq1[i1]);
                    newTrace.alignment2.push_back(seq2[i2]);
                    newTrace.i1--;
                    newTrace.i2--;
                    backtraces.push_back(newTrace);
                }
                backtraces[i].i2--;
            }
            else if ((alignDir & upAlign) == upAlign) {
                backtraces[i].alignment1.push_back(seq1[i1]);
                backtraces[i].alignment2.push_back('-');
                if ((alignDir & diagAlign) == diagAlign) {
                    AlignmentBacktrace newTrace(backtraces[i]);
                    newTrace.alignment1.push_back(seq1[i1]);
                    newTrace.alignment2.push_back(seq2[i2]);
                    newTrace.i1--;
                    newTrace.i2--;
                    backtraces.push_back(newTrace);
                }
                backtraces[i].i1--;
            }
            else { //alignDir contains diagAlign
                backtraces[i].alignment1.push_back(seq1[i1]);
                backtraces[i].alignment2.push_back(seq2[i2]);
                backtraces[i].i1--;
                backtraces[i].i2--;
            }
        }
    }

    double best_score = grid[len1][len2];

    for (int i = 0; i < len1; i++) {
        delete[] grid[i];
        delete[] dir[i];
    }
    delete[] grid[len1];
    delete[] grid;
    delete[] dir;

    std::vector<Alignment> alignments;
    for (int i = 0; i < backtraces.size(); i++) {
        Alignment alignment;
        alignment.score = best_score;
        alignment.sequence1 = backtraces[i].alignment1;
        alignment.sequence2 = backtraces[i].alignment2;
        std::reverse(alignment.sequence1.begin(), alignment.sequence1.end());
        std::reverse(alignment.sequence2.begin(), alignment.sequence2.end());
        alignments.push_back(alignment);
    }
    return alignments;
}
