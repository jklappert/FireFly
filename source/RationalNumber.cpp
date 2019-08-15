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
