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

std::vector<Alignment> GlobalAlignment::backtrace_alignments() {
    unsigned len1 = seq1.size();
    unsigned len2 = seq2.size();

    // for two empty sequences return empty list (as in Bio.pairwise2)
    if (len1 == 0 && len2 == 0) {
        return std::vector<Alignment>();
    }

    /* Tracking back best alignments */
    std::queue<AlignmentBacktrace> backtraces;
    std::vector<Alignment> alignments;

    AlignmentBacktrace initBacktrace;
    initBacktrace.alignment.sequence1 = "";
    initBacktrace.alignment.sequence2 = "";
    initBacktrace.alignment.begin = 0;
    initBacktrace.alignment.end = len1;
    initBacktrace.alignment.score = get_score();
    initBacktrace.i1 = len1-1;
    initBacktrace.i2 = len2-1;
    initBacktrace.trace = trace[len1-1][len2-1];
    DEBUG(std::cout<<"Backtracing... "<<std::endl;)

    backtraces.push(initBacktrace);

    DEBUG(std::cout<<"Backtracing... "<<std::endl;)
    while (!backtraces.empty()) {
        TRACE(std::cout<<"(b:"<<backtraces.size()<<") ";)
        TRACE(std::cout<<"(f:"<<alignments.size()<<") ";)
        auto backtrace = backtraces.front();
        TRACE(std::cout<<"i1,i2: ("<<backtrace.i1<<","<<backtrace.i2<<") ";)
        while (true) {
            TRACE(std::cout<<backtrace.trace<<" ";)
            if ((backtrace.trace & LEFT_ALIGN) == LEFT_ALIGN) {
                backtrace.trace -= LEFT_ALIGN;
                if (backtrace.trace > 0) {
                    AlignmentBacktrace newTrace(backtrace);
                    TRACE(std::cout<<"PUSHING TRACE ("<<newTrace.i1<<","<<newTrace.i2<<") ";)
                    backtraces.push(newTrace);
                }
                backtrace.alignment.sequence1.push_back('-');
                backtrace.alignment.sequence2.push_back(seq2[backtrace.i2]);
                backtrace.i2--;
            }
            else if ((backtrace.trace & UP_ALIGN) == UP_ALIGN) {
                backtrace.trace -= UP_ALIGN;
                if (backtrace.trace > 0) {
                    AlignmentBacktrace newTrace(backtrace);
                    TRACE(std::cout<<"PUSHING TRACE ("<<newTrace.i1<<","<<newTrace.i2<<") ";)
                    backtraces.push(newTrace);
                }
                backtrace.alignment.sequence1.push_back(seq1[backtrace.i1]);
                backtrace.alignment.sequence2.push_back('-');
                backtrace.i1--;
            }
            else { //backtrace.trace contains DIAG_ALIGN
                backtrace.trace -= DIAG_ALIGN;
                if (backtrace.trace > 0) {
                    AlignmentBacktrace newTrace(backtrace);
                    TRACE(std::cout<<"PUSHING TRACE ("<<newTrace.i1<<","<<newTrace.i2<<") ";)
                    backtraces.push(newTrace);
                }
                backtrace.alignment.sequence1.push_back(seq1[backtrace.i1]);
                backtrace.alignment.sequence2.push_back(seq2[backtrace.i2]);
                backtrace.i1--;
                backtrace.i2--;
            }
            if (backtrace.i1 < 0 || backtrace.i2 < 0) {
                break;
            }
            backtrace.trace = trace[backtrace.i1][backtrace.i2];
        }
        DEBUG(std::cout<<"PUSHING ALIGNMENT"<<std::endl;)
        std::reverse(backtrace.alignment.sequence1.begin(), backtrace.alignment.sequence1.end());
        std::reverse(backtrace.alignment.sequence2.begin(), backtrace.alignment.sequence2.end());
        alignments.push_back(backtrace.alignment);
        if (alignments.size() >= MAX_ALIGNMENTS) {
            DEBUG(std::cout<<"REACHED MAX ALIGNMENTS COUNT"<<std::endl;)
            break;
        }
        backtraces.pop();
    }

    return alignments;
}
