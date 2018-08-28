#include "RationalNumber.hpp"

namespace firefly {

  RationalNumber::RationalNumber(mpz_class numerator_, mpz_class denominator_) {
    mpz_class gcd_(gcd(numerator_, denominator_));
    numerator = numerator_ / gcd_;
    denominator = denominator_ / gcd_;
  }

  RationalNumber::RationalNumber() {}

  RationalNumber RationalNumber::operator*(const RationalNumber& rn) {
    mpz_class num(numerator * rn.numerator);
    mpz_class den(denominator * rn.denominator);

    if (den < 0) {
      num = -num;
      den = -den;
    }

    return RationalNumber(num, den);
  }

  std::ostream& operator<< (std::ostream& out, const RationalNumber& a) {
    if (a.denominator == 1) {
      if (a.numerator < 1) {
        out << "(" << a.numerator.get_str() << ")";
      } else {
        out << a.numerator.get_str();
      }
    } else {
      out << "(" << a.numerator.get_str() << "/" << a.denominator.get_str() << ")";
    }

    return out;
  }

}


