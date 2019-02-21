#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "UintHasher.hpp"
#include "Monomial.hpp"
#include "PolynomialFF.hpp"

namespace firefly {
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;

  /**
   * @class Polynomial
   * @brief A container class representing polynomials
   */
  class Polynomial {
  public:
    /**
     *    A constructor for a polynomial with RationalNumber objects as
     *    coefficients
     *    @param coef an unordered map of coefficients and degrees
     */
    Polynomial(const rn_map& coef);
    /**
     *    A constructor for a polynomial by a Monomial object
     *    @param coef a Monomial object
     */
    Polynomial(const Monomial& coef);
    Polynomial();
    Polynomial& operator*=(const RationalNumber&);
    /**
     *  Converts the Polynomial object to a PolynomialFF object
     */
    PolynomialFF convert_to_PolynomialFF();
    /**
     *  Sorts the degress of the polynomial in lexographically order
     */
    void sort();
    /**
     *  Clears the polynomial
     */
    void clear();
    /**
     *  Transforms the Polynomial object to a string where each variable
     *  is replaced by the corresponding symbol in a given vector
     *  @param symbols a vector of symbols, e.g. {"x","y","z"}.
     */
    std::string to_string(const std::vector<std::string>& symbols) const;

    std::vector<Monomial> coefs;
  private:
    uint32_t n;
  };

  std::ostream& operator<< (std::ostream& out, const Polynomial& pol);

}
