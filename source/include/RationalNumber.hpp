#pragma once
#include <gmpxx.h>
#include <iostream>
#include <string>

namespace firefly {

  class RationalNumber {
  public:
    RationalNumber(mpz_class numerator_, mpz_class denominator_);
    RationalNumber();
    RationalNumber operator*(const RationalNumber&);
    RationalNumber& operator+=(const RationalNumber& rn);
    RationalNumber& operator*=(const RationalNumber& rn);
    bool operator==(const RationalNumber&) const;
    RationalNumber operator-();
    std::string string() const;

    mpz_class numerator;
    mpz_class denominator;
  };

  std::ostream& operator<< (std::ostream& out, const RationalNumber&);

}
