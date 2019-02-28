// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================
#include "RationalNumber.hpp"

namespace firefly {

  RationalNumber::RationalNumber(const mpz_class& numerator_, const mpz_class& denominator_) {
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

  RationalNumber& RationalNumber::operator-=(const RationalNumber& rn) {
    if (rn.denominator != denominator) {
      numerator = numerator * rn.denominator - rn.numerator * denominator;
      denominator = denominator * rn.denominator;
    } else numerator -= rn.numerator;

    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }

    mpz_class gcd_(gcd(numerator, denominator));
    numerator = numerator / gcd_;
    denominator = denominator / gcd_;
    return *this;
  }

  RationalNumber& RationalNumber::operator+=(const RationalNumber& rn) {
    if (rn.denominator != denominator) {
      numerator = numerator * rn.denominator + rn.numerator * denominator;
      denominator = denominator * rn.denominator;
    } else numerator += rn.numerator;

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

  RationalNumber RationalNumber::operator-() const {
    return RationalNumber(-numerator, denominator);
  }

  bool RationalNumber::operator==(const RationalNumber& b) const {
    return (numerator == b.numerator && denominator == b.denominator);
  }

  std::string RationalNumber::string() const {
    std::string str;

    if (denominator == 1) {
      if (numerator < 1) {
        str += "(" + numerator.get_str() + ")";
      } else {
        str += numerator.get_str();
      }
    } else {
      if (numerator < 1) {
        str += "(" + numerator.get_str() + "/" + denominator.get_str() + ")";
      } else {
        str += numerator.get_str() + "/" + denominator.get_str();
      }
    }

    return str;
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
