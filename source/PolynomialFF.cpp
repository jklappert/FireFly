#include "PolynomialFF.hpp"

namespace firefly {

  PolynomialFF::PolynomialFF(std::vector<FFInt> coef_) : coef(coef_) {
    deg = coef.size() - 1;
  }

  PolynomialFF::PolynomialFF() {}

  FFInt PolynomialFF::calc(FFInt x) {
    const uint64_t prime = coef.at(0).p;
    FFInt res(0, prime);

    for (uint i = 0; i < (uint) coef.size(); i++) {
      FFInt exp(i, prime);
      res += coef.at(i) * x.pow(exp);
    }

    return res;
  }

  PolynomialFF PolynomialFF::operator+(const PolynomialFF &b) {
    PolynomialFF a = *this;
    std::vector<FFInt> newCoefs {};

    if (a.deg >= b.deg) {
      newCoefs = a.coef;

      for (int i = 0; i <= b.deg; i++) newCoefs.at(i) += b.coef.at(i);
    } else {
      newCoefs = b.coef;

      for (int i = 0; i <= a.deg; i++) newCoefs.at(i) += a.coef.at(i);
    }

    return PolynomialFF(newCoefs);
  }

  PolynomialFF PolynomialFF::operator-(const PolynomialFF &b) {
    PolynomialFF a = *this;
    std::vector<FFInt> newCoefs {};

    if (a.deg >= b.deg) {
      newCoefs = a.coef;

      for (int i = 0; i <= b.deg; i++) newCoefs.at(i) -= b.coef.at(i);
    } else {
      newCoefs = b.coef;
      FFInt zero(0, a.coef.at(1).p);

      for (int i = 0; i <= a.deg; i++) newCoefs.at(i) = a.coef.at(i) - newCoefs.at(i);

      for (int i = a.deg + 1; i <= b.deg; i++) newCoefs.at(i) = zero - newCoefs.at(i);
    }

    return PolynomialFF(newCoefs);
  }

  PolynomialFF PolynomialFF::operator*(const PolynomialFF &b) {
    PolynomialFF a = *this;
    std::vector<FFInt> newCoefs {};
    const double newDeg = a.deg + b.deg;

    for (int i = 0; i <= newDeg; i++) {
      FFInt ffint(0, a.coef.at(0).p);

      for (int j = 0; j <= i; j++) {
        if (a.deg >= j && b.deg >= (i - j)) {
          ffint += a.coef.at(j) * b.coef.at(i - j);
        }
      }

      newCoefs.push_back(ffint);
    }

    return PolynomialFF(newCoefs);
  }

  PolynomialFF &PolynomialFF::operator=(const PolynomialFF &a) {
    coef = a.coef;
    deg = a.deg;
    return *this;
  }

  PolynomialFF PolynomialFF::operator*(const FFInt &a) {
    std::vector<FFInt> newCoefs {};

    for (auto coefficient : coef) {
      newCoefs.push_back(coefficient * a);
    }

    return PolynomialFF(newCoefs);
  }

  PolynomialFF PolynomialFF::operator/(const FFInt &a) {
    std::vector<FFInt> newCoefs {};

    for (auto coefficient : coef) {
      newCoefs.push_back(coefficient / a);
    }

    return PolynomialFF(newCoefs);
  }



  std::ostream &operator<<(std::ostream &out, const PolynomialFF &a) {
    if (a.coef.size() == 1) return out << a.coef.at(0).n;

    for (int i = 0; i < (int) a.coef.size(); i++) {
      const uint64_t n = a.coef.at(i).n;

      if (i == 0) {
        out << n << "*x^" << i;
      } else if (n != 0) {
        out << " + " << n << "*x^" << i;
      }
    }

    return out;
  }

}