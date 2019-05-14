//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

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
