#ifndef MATCH_SCORERS_hpp
#define MATCH_SCORERS_hpp

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
        virtual double getScore(char el1, char el2) const = 0;
};

/*! A match score is the score of identical elements, otherwise mismatch score */
class ConstMatchScorer : public MatchScorer {
    private:
        double matchScore;
        double mismatchScore;
    public:
        ConstMatchScorer(double, double);
        double getScore(char el1, char el2) const;
};

/*! A match score is based on dictionary, where each pair of characters has its score */
class DictMatchScorer : public MatchScorer {
    private:
        std::map<char, int> elementsIndices;
        double** scoreMatrix;
    public:
        DictMatchScorer(std::map<char, int>, double**);
        double getScore(char el1, char el2) const;
};

/*! A match score is a result of a callback function */
class CallbackMatchScorer : public MatchScorer {
    private:
        double(*scoreCallback)(char, char);
    public:
        CallbackMatchScorer(double(*)(char, char));
        double getScore(char el1, char el2) const;
};


#endif
