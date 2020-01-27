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
#include "HornerGenerator.hpp"

namespace firefly {

  RationalFunction::RationalFunction(const Polynomial& n, const Polynomial& d) : numerator(n), denominator(d) {}

  RationalFunction::RationalFunction() {}

  std::string RationalFunction::to_string(const std::vector<std::string>& vars) const {
    std::string str = "";
    if (!factors.empty()) {
      std::string num_str = "(";
      std::string den_str = "(";

      for (const auto& factor : factors) {
        if (!factor.numerator.coefs.empty())
          num_str += "(" + factor.numerator.to_string(vars) + ")*";

        if (!factor.denominator.coefs.empty())
          den_str += "(" + factor.denominator.to_string(vars) + ")*";
      }

      num_str += "(" + numerator.to_string(vars) + "))/";

      if(!denominator.coefs.empty())
        den_str += "(" + denominator.to_string(vars) + "))";
      else
        den_str += "(1))";

      str += num_str + den_str;
    } else {
      str += "(" + numerator.to_string(vars) + ")/(";
      if(!denominator.coefs.empty())
        str += denominator.to_string(vars) + ")";
      else
        str += "1)";
    }

    return str;
  }

  std::ostream& operator<< (std::ostream& out, const RationalFunction& rf) {
    out << "Numerator: " << rf.numerator;
    out << "Denominator: " << rf.denominator;
    return out;
  }

  std::string RationalFunction::generate_horner(const std::vector<std::string>& vars) const {
    return "(" + numerator.generate_horner(vars) + ")/(" + denominator.generate_horner(vars) + ")";
  }

  void RationalFunction::add_factor(const RationalFunction& factor) {
    factors.reserve(factors.size() + 1);
    factors.emplace_back(factor);
  }

  std::vector<RationalFunction> RationalFunction::get_factors() const {
    return factors;
  }
}
