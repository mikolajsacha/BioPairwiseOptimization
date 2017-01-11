#include "../include/global_alignment.h"
#include "../include/match_scorers.h"
#include "../include/model.h"
#include "../include/debug.h"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <iostream>

GlobalAlignment::GlobalAlignment(const std::string& seq1, const std::string& seq2) : seq1(seq1), seq2(seq2) {
    alignments_count = 0;
    unsigned len1 = seq1.size();
    unsigned len2 = seq2.size();

    grid = new float*[len1 + 1];
    trace = new char*[len1];
    for (unsigned i = 0; i < len1; i++) {
        grid[i] = new float[len2 + 1];
        trace[i] = new char[len2];
    }
    grid[len1] = new float[len2 + 1];
}

GlobalAlignment::~GlobalAlignment() {
    unsigned len1 = seq1.size();
    for (unsigned i = 0; i < len1; i++) {
        delete[] grid[i];
        delete[] trace[i];
    }
    delete[] grid[len1];
    delete[] grid;
    delete[] trace;
}

float GlobalAlignment::get_score() {
    return grid[seq1.size()][seq2.size()];
}

void GlobalAlignment::populate_matrix_linear_gap_penalty_only_grid(MatchScorer* scorer, float penalty) {
    unsigned len1 = seq1.size();
    unsigned len2 = seq2.size();

    for (unsigned i = 0; i <= len1; i++) {
        grid[i][0] = i * penalty;
    }
    for (unsigned i = 0; i <= len2; i++) {
        grid[0][i] = i * penalty;
    }

    DEBUG(std::cout<<"Calculating only grid matrix..."<<std::endl;)
    for (unsigned i = 1; i <= len1; i++) {
        for (unsigned j = 1; j <= len2; j++) {
            float matchScore = scorer->getScore(seq1[i-1], seq2[j-1]);
            float diagScore = grid[i-1][j-1] + matchScore;
            float upScore = grid[i-1][j] + penalty;
            float leftScore = grid[i][j-1] + penalty;
            TRACE(std::cout<<"match: "<<matchScore<<", diag: "<<diagScore<<;)
            TRACE(std::cout<<", up: "<<upScore<<", left: "<<leftScore<<std::endl;)

            if (diagScore >= upScore) {
                if (diagScore >= leftScore) { 
                    grid[i][j] = diagScore;
                }
                else {
                    grid[i][j] = leftScore;
                }
            }
            else { //diagScore < upScore
                if (upScore >= leftScore) {
                    grid[i][j] = upScore;
                }
                else {
                    grid[i][j] = leftScore;
                }
            }
        }
    }
}

void GlobalAlignment::populate_matrix_linear_gap_penalty(MatchScorer* scorer, float penalty) {
    unsigned len1 = seq1.size();
    unsigned len2 = seq2.size();

    for (unsigned i = 0; i <= len1; i++) {
        grid[i][0] = i * penalty;
    }
    for (unsigned i = 0; i <= len2; i++) {
        grid[0][i] = i * penalty;
    }

    DEBUG(std::cout<<"Calculating grid and trace matrices..."<<std::endl;)
    for (unsigned i = 1; i <= len1; i++) {
        for (unsigned j = 1; j <= len2; j++) {
            float matchScore = scorer->getScore(seq1[i-1], seq2[j-1]);
            float diagScore = grid[i-1][j-1] + matchScore;
            float upScore = grid[i-1][j] + penalty;
            float leftScore = grid[i][j-1] + penalty;
            TRACE(std::cout<<"match: "<<matchScore<<", diag: "<<diagScore<<;)
            TRACE(std::cout<<", up: "<<upScore<<", left: "<<leftScore<<std::endl;)

            if (diagScore > upScore) {
                if (diagScore > leftScore) { //upScore < leftScore < diagScore
                    grid[i][j] = diagScore;
                    trace[i-1][j-1] = DIAG_ALIGN;
                }
                else if (diagScore == leftScore) { // upScore < diagScore == leftScore
                    grid[i][j] = diagScore;
                    trace[i-1][j-1] = DIAG_LEFT_ALIGN;
                }
                else { // upScore < diagScore < leftScore
                    grid[i][j] = leftScore;
                    trace[i-1][j-1] = LEFT_ALIGN;
                }
            }
            else if (diagScore == upScore) {
                if (diagScore > leftScore) { //leftScore < upScore == diagScore
                    grid[i][j] = diagScore;
                    trace[i-1][j-1] = DIAG_UP_ALIGN;
                }
                else if (diagScore == leftScore) { // upScore == diagScore == leftScore
                    grid[i][j] = diagScore;
                    trace[i-1][j-1] = ALL_ALIGN;
                }
                else { // upScore == diagScore < leftScore
                    grid[i][j] = leftScore;
                    trace[i-1][j-1] = LEFT_ALIGN;
                }
            }
            else { // diagScore < upScore
                if (upScore > leftScore) { // diagScore < upScore, leftScore < upScore
                    grid[i][j] = upScore;
                    trace[i-1][j-1] = UP_ALIGN;
                }
                else if (upScore == leftScore) { // diagScore < upScore == leftScore
                    grid[i][j] = upScore;
                    trace[i-1][j-1] = UP_LEFT_ALIGN;
                }
                else { // diagScore < upScore < leftScore
                    grid[i][j] = leftScore;
                    trace[i-1][j-1] = LEFT_ALIGN;
                }
            }
        }
    }
}

void GlobalAlignment::check_for_new_backtrace(unsigned i) {
    if (backtraces[i].trace > 0 && alignments_count < MAX_ALIGNMENTS) {
        alignments[alignments_count].sequence1 = alignments[i].sequence1;
        alignments[alignments_count].sequence2 = alignments[i].sequence2;
        alignments[alignments_count].score = alignments[i].score;
        alignments[alignments_count].begin = alignments[i].begin;
        alignments[alignments_count].end = alignments[i].end;

        backtraces[alignments_count].i1 = backtraces[i].i1;
        backtraces[alignments_count].i2 = backtraces[i].i2;
        backtraces[alignments_count].trace = backtraces[i].trace;
        alignments_count++;
        TRACE(std::cout<<"PUSHING TRACE ("<<backtraces[i].i1<<","<<backtraces[i].i2<<") ";)
    }
}

void GlobalAlignment::backtrace_alignments() {
    unsigned len1 = seq1.size();
    unsigned len2 = seq2.size();

    // for two empty sequences return empty list (as in Bio.pairwise2)
    if (len1 == 0 && len2 == 0) {
        alignments_count = 0;
        return;
    }

    /* Tracking back best alignments */
    Alignment initAlignment;
    initAlignment.sequence1 = "";
    initAlignment.sequence2 = "";
    initAlignment.begin = 0;
    initAlignment.end = len1;
    initAlignment.score = get_score();

    AlignmentBacktrace initBacktrace;
    initBacktrace.i1 = len1-1;
    initBacktrace.i2 = len2-1;
    initBacktrace.trace = trace[len1-1][len2-1];

    backtraces[0] = initBacktrace;
    alignments[0] = initAlignment;
    alignments_count = 1;

    DEBUG(std::cout<<"Backtracing... "<<std::endl;)
    for (unsigned i= 0; i < alignments_count; i++) {
        TRACE(std::cout<<"(a:"<<alignments_count<<") ";)
        TRACE(std::cout<<"i1,i2: ("<<alignments[i].i1<<","<<alignments[i].i2<<") ";)
        while (true) {
            TRACE(std::cout<<alignments[i].trace<<" ";)
            if ((backtraces[i].trace & LEFT_ALIGN) == LEFT_ALIGN) {
                backtraces[i].trace -= LEFT_ALIGN;
                check_for_new_backtrace(i);
                alignments[i].sequence1.push_back('-');
                alignments[i].sequence2.push_back(seq2[backtraces[i].i2]);
                backtraces[i].i2--;
            }
            else if ((backtraces[i].trace & UP_ALIGN) == UP_ALIGN) {
                backtraces[i].trace -= UP_ALIGN;
                check_for_new_backtrace(i);
                alignments[i].sequence1.push_back(seq1[backtraces[i].i1]);
                alignments[i].sequence2.push_back('-');
                backtraces[i].i1--;
            }
            else { //backtraces[i].trace contains DIAG_ALIGN
                backtraces[i].trace -= DIAG_ALIGN;
                check_for_new_backtrace(i);
                alignments[i].sequence1.push_back(seq1[backtraces[i].i1]);
                alignments[i].sequence2.push_back(seq2[backtraces[i].i2]);
                backtraces[i].i1--;
                backtraces[i].i2--;
            }
            if (backtraces[i].i1 < 0 || backtraces[i].i2 < 0) {
                break;
            }
            backtraces[i].trace = trace[backtraces[i].i1][backtraces[i].i2];
        }
        DEBUG(std::cout<<"FINISHED BACKTRACING ALIGNMENT"<<std::endl;)
        std::reverse(alignments[i].sequence1.begin(), alignments[i].sequence1.end());
        std::reverse(alignments[i].sequence2.begin(), alignments[i].sequence2.end());
    }
    DEBUG(
        if (alignments_count >= MAX_ALIGNMENTS) {
            std::cout<<"REACHED MAX ALIGNMENTS COUNT"<<std::endl
        }
    )
}
