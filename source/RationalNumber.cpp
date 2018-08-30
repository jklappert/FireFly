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

    mpz_class gcd_(gcd(num, den));
    numerator = num / gcd_;
    denominator = den / gcd_;
    return *this;
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

  RationalNumber& RationalNumber::operator+=(const RationalNumber& rn) {
    if (rn.denominator != denominator) {
      numerator = numerator * rn.denominator + rn.numerator * denominator;
      denominator = denominator * rn.denominator;
    }
    else numerator += rn.numerator;

    mpz_class gcd_(gcd(numerator, denominator));
    numerator = numerator / gcd_;
    denominator = denominator / gcd_;
    return *this;
  }

  RationalNumber& RationalNumber::operator*=(const RationalNumber& rn) {
    mpz_class num(numerator * rn.numerator);
    mpz_class den(denominator * rn.denominator);

    if (den < 0) {
      num = -num;
      den = -den;
    }

    mpz_class gcd_(gcd(num, den));
    numerator = num / gcd_;
    denominator = den / gcd_;
    return *this;
  }

  RationalNumber RationalNumber::operator-() {
    numerator = - numerator;
    return *this;
  }

  bool RationalNumber::operator==(const RationalNumber& b) const {
    return (numerator == b.numerator && denominator == b.denominator);
  }
}


