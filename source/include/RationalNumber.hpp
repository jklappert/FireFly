#pragma once
#include <gmpxx.h>
#include <iostream>

namespace firefly {
  
  class RationalNumber {
  public:
    RationalNumber(mpz_class nominator_, mpz_class denominator_);
    
    mpz_class nominator;
    mpz_class denominator;
  };
  
  std::ostream &operator<< (std::ostream &out, const RationalNumber &a);
  
}