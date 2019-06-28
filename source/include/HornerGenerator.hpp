#pragma once

#include "Polynomial.hpp"

namespace firefly {

  /**
   *  Generates a Horner form of a polynomial over a finite field
   *  @param monomials A ff_map representing a polynomial over a finite field
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_ff(const ff_map& monomials, const std::vector<std::string>& vars, uint32_t index = 0);
  /**
   *  Generates a Horner form of a polynomial over the rationals
   *  @param monomials A rn_map representing a polynomial over the rationals
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_rn(const rn_map& monomials, const std::vector<std::string>& vars, uint32_t index = 0);
  /**
   *  Generates a Horner form of a polynomial over the rationals
   *  @param monomials A vector of monomial objects representing a polynomial over the rationals
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_mon(const std::vector<Monomial> monomials, const std::vector<std::string>& vars, uint32_t index = 0);
}
