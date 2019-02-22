// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

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
    RationalNumber(mpz_class numerator_, mpz_class denominator_);
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
