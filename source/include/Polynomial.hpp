#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "UintHasher.hpp"
#include "Monomial.hpp"

namespace firefly {
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;

  class Polynomial {
  public:
    /**
     *    A constructor for a polynomial with RationalNumber objects as
     *    coefficients
     *    @param coefs_ an unordered map of coefficients and degrees
     */
    Polynomial(const rn_map& coef);
    Polynomial();
    Polynomial operator*(const RationalNumber&);
    Polynomial operator*(const Monomial&);
    Polynomial operator*(const Polynomial&);
    Polynomial& operator+=(const Polynomial&);
    Polynomial& operator+=(const Monomial&);
    Polynomial homogenize(uint degree);
    void sort();
    void clear();
    std::string string(const std::vector<std::string>& symbols) const;

    std::vector<Monomial> coefs;  /**< The vector which holds all coefficients*/
  };

  std::ostream& operator<< (std::ostream& out, const Polynomial& pol);

}
