#ifndef MATCH_SCORERS_H
#define MATCH_SCORERS_H

#include <map>

/*! Abstract class containing method which returns alignment score for two sequence elements */
class MatchScorer {
    public:
        //! Returns score for two characters (sequences' elements)
        /*!
          \param el1 an element from the first sequence
          \param el2 an element from the second sequence
          \return alignment score for given elements
        */
        virtual float getScore(char el1, char el2) const = 0;
        virtual ~MatchScorer() {}
};

/*! A match score is the score of identical elements, otherwise mismatch score */
class ConstMatchScorer : public MatchScorer {
    private:
        float matchScore;
        float mismatchScore;
    public:
        ConstMatchScorer(float, float);
        float getScore(char el1, char el2) const;
        ConstMatchScorer() {}
};

/*! A match score is based on dictionary, where each pair of characters has its score */
class DictMatchScorer : public MatchScorer {
    private:
        std::map<char, int> elementsIndices;
        float** scoreMatrix;
    public:
        DictMatchScorer(std::map<char, int>, float**);
        float getScore(char el1, char el2) const;
        DictMatchScorer() {}
};

/*! A match score is a result of a callback function */
class CallbackMatchScorer : public MatchScorer {
    private:
        float(*scoreCallback)(char, char);
    public:
        CallbackMatchScorer(float(*)(char, char));
        float getScore(char el1, char el2) const;
        CallbackMatchScorer() {}
};


#endif
