#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"
#include "Polynomial.hpp"

namespace firefly {

  class PolyReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters as an integer
     */
    PolyReconst (int n_);
    /**
     *    Calls the reconstruction algorithm
     *    @returns a vector with a pair of integers corrisponding to the
     *    coefficient of the polynom. The vector is ordered in an anscending
     *    way such that the first coefficient corresponds to to z^0,...
     */
    std::vector<FFInt> reconst();
    void constrCanonical();
    Polynomial canonical;
  private:
    /**
     *    A numerical black box function which provides the reconstruction
     *    algorithm with the finite field member f(y)
     *    @param p the prime number which defines the finite field
     *    @param y a member of the finite field which should be evaluated in f(y)
     */
    FFInt num (uint64_t p, const FFInt &y);
    /**
     *    Computes the coefficient a_i recursevly using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param i the order of a_i
     *    @param ip recursion order
     *    @param num f(y_i)
     *    @returns a_i
     */
    FFInt compAi (int i, int ip, const FFInt &num);
    Polynomial iterateCanonical (uint i);

    int n; /**< The number of parameters */
    const uint64_t prime; /**< The prime number defining the finite field */
    std::vector<FFInt> ai {}; /**< A vector which holds all coefficients a_i */
    std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
  };
}
