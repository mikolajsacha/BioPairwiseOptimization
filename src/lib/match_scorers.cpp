#include "../include/match_scorers.h"
#include <iostream>

ConstMatchScorer::ConstMatchScorer(float matchScore, float mismatchScore)
    : matchScore(matchScore), mismatchScore(mismatchScore) {}

float ConstMatchScorer::getScore(char el1, char el2) const {
    return el1 == el2 ? this->matchScore : this->mismatchScore;
}

DictMatchScorer::DictMatchScorer(std::map<char, int> elementsIndices, float** scoreMatrix)
    : elementsIndices(elementsIndices), scoreMatrix(scoreMatrix) {}

float DictMatchScorer::getScore(char el1, char el2) const {
    int i1 = this->elementsIndices.at(el1);
    int i2 = this->elementsIndices.at(el2);
    return this->scoreMatrix[i1][i2];
}

CallbackMatchScorer::CallbackMatchScorer(float(*scoreCallback)(char, char))
    : scoreCallback(scoreCallback) {}

float CallbackMatchScorer::getScore(char el1, char el2) const {
    return this->scoreCallback(el1, el2);
}
