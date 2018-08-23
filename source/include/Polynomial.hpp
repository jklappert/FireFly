#pragma once

#include <vector>
#include <unordered_map>
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
    std::vector<Monomial> coefs;  /**< The vector which holds all coefficients*/
    Polynomial operator*(const RationalNumber &);
  };

  std::ostream &operator<< (std::ostream &out, const Polynomial &pol);

}
