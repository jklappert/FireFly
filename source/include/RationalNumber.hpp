#pragma once
#include <gmpxx.h>
#include <iostream>

namespace firefly {

  class RationalNumber {
  public:
    RationalNumber(mpz_class numerator_, mpz_class denominator_);
    RationalNumber();
    mpz_class numerator;
    mpz_class denominator;
    RationalNumber operator*(const RationalNumber &);
  };

  std::ostream &operator<< (std::ostream &out, const RationalNumber &a);

}
