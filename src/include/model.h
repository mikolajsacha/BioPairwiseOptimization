#ifndef MODEL_H
#define MODEL_H

#include <string>

/*! An Alignment result. Matches format of results in Bio.pairwise2 module */
struct Alignment {
    std::string sequence1; /*!< First sequence in aligned format (with hyphens as gaps) */
    std::string sequence2; /*!< Second sequence in aligned format (with hyphens as gaps) */
    double score; /*! Score for this alignment */
    int begin; /*!< Index where alignment begins */
    int end; /*!< Index where alignment ends */
};

/*! Enum representing direction of possible best alignment, when backtracing alignment in Needleman-Wunsh algorithm */
enum AlignmentDirection {
    leftAlign = 1,     // 0b00000001
    upAlign = 2,       // 0b00000010
    upLeftAlign = 3,   // 0b00000011
    diagAlign = 4,     // 0b00000100
    diagLeftAlign = 5, // 0b00000101 contains left and diagonal
    diagUpAlign = 6,   // 0b00000110 contains left and up
    allAlign = 7       // 0b00000111 contains all directions
};

/*! A structure representing state of tracing back an alignment */
struct AlignmentBacktrace {
    std::string alignment1; /*!< Current alignment of first sequence */
    std::string alignment2; /*!< Current alignment of second sequence */
    int i1; /*!< Current index on first sequence axis */
    int i2; /*!< Current index on second sequence axis */
    bool isFinished; /*!< Marks that current alignment has finished backtracing */
};

#endif
