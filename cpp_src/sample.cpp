#include <iostream>
#include "alignment.h"
#include "match_scorers.h"

int main() {
    ConstMatchScorer a(1.0, -1.0);
    std::cout<<a.getScore('a', 'b')<<endl;
}
