#include "RationalFunction.hpp"

namespace firefly {

  RationalFunction::RationalFunction(Polynomial n, Polynomial d) : numerator(n), denominator(d) {}

  std::ostream &operator<< (std::ostream &out, const RationalFunction &rf) {
    out << "Numerator: " << rf.numerator;
    out << "Denominator: " << rf.denominator;
    return out;
  }
}
