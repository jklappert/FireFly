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

#pragma once

#include <gmpxx.h>
#include <iostream>
#include <string>

namespace firefly {
  /**
   * @class RationalNumber
   * @brief A container class representing rational numbers
   */
  class RationalNumber {
  public:
    /**
     *  Constructor of a RationalNumber object
     *  @param numberator_ the numerator as a mpz_class
     *  @param denominator_ the denominator as a mpz_class
     */
    RationalNumber(const mpz_class& numerator_, const mpz_class& denominator_);
    RationalNumber();
    RationalNumber operator*(const RationalNumber&);
    RationalNumber& operator+=(const RationalNumber& rn);
    RationalNumber& operator-=(const RationalNumber& rn);
    RationalNumber& operator*=(const RationalNumber& rn);
    bool operator==(const RationalNumber&) const;
    RationalNumber operator-() const;
    std::string string() const;

    mpz_class numerator;
    mpz_class denominator;
  };

  std::ostream& operator<< (std::ostream& out, const RationalNumber&);

}
