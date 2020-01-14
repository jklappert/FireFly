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

#include "RationalFunctionFF.hpp"

namespace firefly {

  RationalFunctionFF::RationalFunctionFF(const PolynomialFF& n, const PolynomialFF& d) : numerator(n), denominator(d) {}

  RationalFunctionFF::RationalFunctionFF() {}

  std::string RationalFunctionFF::to_string(const std::vector<std::string>& vars) const {
    std::string str = "(" + numerator.to_string(vars) + ")/(";
    if(!denominator.coefs.empty())
      str += denominator.to_string(vars) + ")";
    else
      str += "1)";
    return str;
  }

  std::ostream& operator<< (std::ostream& out, const RationalFunctionFF& rf) {
    out << "Numerator: " << rf.numerator;
    out << "Denominator: " << rf.denominator;
    return out;
  }
}
