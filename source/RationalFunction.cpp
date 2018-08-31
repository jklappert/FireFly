#include "RationalFunction.hpp"

namespace firefly {

  RationalFunction::RationalFunction(Polynomial n, Polynomial d) : numerator(n), denominator(d) {}

  RationalFunction::RationalFunction() {}

  std::string RationalFunction::string(const std::vector<std::string>& symbols) const {
    std::string str = "(" + numerator.string(symbols) + ")/(" + denominator.string(symbols) + ")";
    return str;
  }

  std::ostream& operator<< (std::ostream& out, const RationalFunction& rf) {
    out << "Numerator: " << rf.numerator;
    out << "Denominator: " << rf.denominator;
    return out;
  }

}
