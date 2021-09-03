//
// Created by Ale on 21/06/2021.
//

#include "../include/modelica/utils/VariableTracker.h"

VariableTracker::VariableTracker(const string &name, const bool isArray, const bool isDerivative, const uint16_t dim)
        : _name(name), _isArray(isArray), _isDerivative(isDerivative), _dim(dim) {

}

void VariableTracker::setRanges(const list<Range> &ranges) {
    _ranges = ranges;
}

const string &VariableTracker::getName() const {
    return _name;
}

const bool VariableTracker::getIsArray() const {
    return _isArray;
}

const bool VariableTracker::getIsDerivative() const {
    return _isDerivative;
}

const uint16_t VariableTracker::getDim() const {
    return _dim;
}

const list<Range> &VariableTracker::getRanges() const {
    return _ranges;
}

void VariableTracker::dump(void) const {
    printf("\n* 💾Variable:\n");
    cout << "\tname: " << getName();
    cout << "\tderivative: " << (getIsDerivative() ? " yes" : " no");
    cout << "\tisArray: " << (getIsArray() ? " yes" : " no") << endl;
    if (getIsArray()) {
        cout << "\t SIZE: " << getDim();
        for (const auto &item : getRanges()) {
            cout << "\n\t\trange [" << item.leftValue << "," << item.rightValue << "];";
        }
        printf("\n");
    }
}


Range::Range(int leftValue, int rightValue) : leftValue(leftValue), rightValue(rightValue) {}