#include "benchmarker/residue_name.hpp"
#include <cstring>

using namespace benchmarker;

static char check_digit_ (char c);
static char clamp_(char c);

ResidueName::ResidueName(const std::string& s) :
    super({{0, 0, 0}}) {

    switch (s.length()) {
        case 3:
            (*this)[2] = check_digit_(s[2]);
        case 2:
            (*this)[1] = check_digit_(s[1]);
        case 1:
            (*this)[0] = check_digit_(s[0]);
            break;
        default:
            throw std::length_error("Cannot have a residue name with given character length");
            break;
    }
}

ResidueName::ResidueName(const char* s) :
    super({{0, 0, 0}}) {

    switch (std::strlen(s)) {
        case 3:
            (*this)[2] = check_digit_(s[2]);
        case 2:
            (*this)[1] = check_digit_(s[1]);
        case 1:
            (*this)[0] = check_digit_(s[0]);
            break;
        default:
            throw std::length_error("Cannot have a residue name with given character length");
            break;
    }
}

char check_digit_(char c) {
#ifndef NDEBUG
    if ( !((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z')) ) {
        throw std::range_error("Invalid character");
    }
#endif
    return c;
}

char clamp_(char c) {

    if (c == 0) {
        return 36;
    }

    if (c >= '0' && c <= '9') {
        return c - '0';
    }
    
    return 10 + c - 'A';
}

size_t ResidueName::hash() const {
    return clamp_((*this)[0]) + clamp_((*this)[1]) * 37 + clamp_((*this)[2]) * (37 * 37);
}

const std::array<char, 3>& ResidueName::operator*() const {
    return *this;
}

size_t ResidueNameHash::operator()(const ResidueName& resn) const {
    return resn.hash();
}
