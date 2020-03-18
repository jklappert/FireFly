//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "Monomial.hpp"

namespace firefly {

  Monomial::Monomial(const std::vector<uint32_t>& powers_, const RationalNumber& coef_) : powers(powers_), coef(coef_) {}

  bool Monomial::operator<(const Monomial& m2) {
    uint32_t deg1 = 0;
    uint32_t deg2 = 0;

    for (uint32_t i = 0; i < powers.size(); i++) {
      deg1 += powers[i];
      deg2 += m2.powers[i];
    }

    if (deg1 < deg2) {
      return true;
    } else {
      // since (0,0,1,0) > (0,1,0,0) we need to change the relational operator
      return deg1 == deg2 && powers > m2.powers;
    }
  }

  bool Monomial::operator>(const Monomial& m2) {
    uint32_t deg1 = 0;
    uint32_t deg2 = 0;

    for (uint32_t i = 0; i < powers.size(); ++i) {
      deg1 += powers[i];
      deg2 += m2.powers[i];
    }

    if (deg1 > deg2) {
      return true;
    } else {
      // since (0,0,1,0) > (0,1,0,0) we need to change the relational operator
      return deg1 == deg2 && powers < m2.powers;
    }
  }

  Monomial Monomial::operator*(const Monomial& b) {
    Monomial a = *this;
    a.coef = a.coef * b.coef;

    for (uint32_t i = 0; i < powers.size(); ++i) {
      a.powers[i] += b.powers[i];
    }

    return a;
  }

  Monomial Monomial::operator-() const {
    return Monomial(powers, - coef);
  }
}
