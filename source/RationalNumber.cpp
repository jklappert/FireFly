#include "RationalNumber.hpp"

namespace firefly {

  RationalNumber::RationalNumber(mpz_class numerator_, mpz_class denominator_) {
    mpz_class gcd_(gcd(numerator_, denominator_));
    numerator = numerator_ / gcd_;
    denominator = denominator_ / gcd_;
  }

  std::ostream &operator<< (std::ostream &out, const RationalNumber &a) {
    out << "(" << a.numerator.get_str() << "/" << a.denominator.get_str() << ")";
    return out;
  }

}
