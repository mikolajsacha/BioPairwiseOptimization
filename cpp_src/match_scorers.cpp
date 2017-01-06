#include "match_scorers.h"
#include <iostream>

ConstMatchScorer::ConstMatchScorer(double matchScore, double mismatchScore)
    : matchScore(matchScore), mismatchScore(mismatchScore) {}

double ConstMatchScorer::getScore(char el1, char el2) const {
    return el1 == el2 ? this->matchScore : this->mismatchScore;
}

DictMatchScorer::DictMatchScorer(std::map<char, int> elementsIndices, double** scoreMatrix)
    : elementsIndices(elementsIndices), scoreMatrix(scoreMatrix) {}

double DictMatchScorer::getScore(char el1, char el2) const {
    int i1 = this->elementsIndices.at(el1);
    int i2 = this->elementsIndices.at(el2);
    return this->scoreMatrix[i1][i2];
}

CallbackMatchScorer::CallbackMatchScorer(double(*scoreCallback)(char, char))
    : scoreCallback(scoreCallback) {}

double CallbackMatchScorer::getScore(char el1, char el2) const {
    return this->scoreCallback(el1, el2);
}
