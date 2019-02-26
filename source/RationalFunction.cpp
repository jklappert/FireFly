// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================
#include "RationalFunction.hpp"

namespace firefly {

  RationalFunction::RationalFunction(const Polynomial& n, const Polynomial& d) : numerator(n), denominator(d) {}

  RationalFunction::RationalFunction() {}

  std::string RationalFunction::to_string(const std::vector<std::string>& symbols) const {
    std::string str = "(" + numerator.to_string(symbols) + ")/(" + denominator.to_string(symbols) + ")";
    return str;
  }

  std::ostream& operator<< (std::ostream& out, const RationalFunction& rf) {
    out << "Numerator: " << rf.numerator;
    out << "Denominator: " << rf.denominator;
    return out;
  }

}
