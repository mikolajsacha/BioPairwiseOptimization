#ifndef MODEL_H
#define MODEL_H

#define MAX_ALIGNMENTS 1000 // same as in bio.pairwise2

/*! Alignment directions (in code they are of type 'char' for space efficiency) */
#define LEFT_ALIGN 1
#define UP_ALIGN 2
#define UP_LEFT_ALIGN 3
#define DIAG_ALIGN 4
#define DIAG_LEFT_ALIGN 5
#define DIAG_UP_ALIGN 6
#define ALL_ALIGN 7

#include <string>

/*! An Alignment result. Matches format of results in Bio.pairwise2 module */
struct Alignment {
    std::string sequence1; /*!< First sequence in aligned format (with hyphens as gaps) */
    std::string sequence2; /*!< Second sequence in aligned format (with hyphens as gaps) */
    float score; /*! Score for this alignment */
    unsigned begin; /*!< Index where alignment begins */
    unsigned end; /*!< Index where alignment ends */
};


/*! A structure representing state of tracing back an alignment */
struct AlignmentBacktrace {
    int i1; /*!< Current index on first sequence axis */
    int i2; /*!< Current index on second sequence axis */
    char trace; /*< Cached value from trace matrix */
};

#endif
