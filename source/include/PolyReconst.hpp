#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>
#include "Polynomial.hpp"
#include "RationalNumber.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint>, mpz_class, UintHasher> mpz_map;
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;

  class PolyReconst {
  public:
    /**
     *    A constructor
     *    @param n_ The number of parameters as an integer
     */
    PolyReconst(uint n_, const std::vector<FFInt> &anchor_points = std::vector<FFInt> (), const int deg_inp = -1);
    /**
     *    Default constructor. Should not be used explicitly.
     */
    PolyReconst();
    /**
     *    Calls the reconstruction algorithm
     *    @return the reconstructed Polynomial
     *    @throw runtimeerror if the prime numbers are not sufficient to reconstruct
     *    rational coefficients
     */
    void feed(const std::vector<FFInt>& yis, const FFInt& num);
    /**
     * 
     */
    void feed(const std::vector<std::vector<uint>>& degs, const std::vector<FFInt>& new_yis, const FFInt& num);
    bool done = false;
    bool new_prime = false;
    uint next_zi = 1;
    uint prime_number = 0;
    Polynomial get_result();
    std::vector<uint> curr_zi_order{};
  private:
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param zi the integer i to a zi
     *    @param ai The vector of previously computed ai
     *    @param num f(y_i)
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @return a(i)
     */
    PolynomialFF comp_ai(const uint zi, int i, int ip, const PolynomialFF& num, std::vector<PolynomialFF>& ai);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param zi the integer i to a zi
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @return The vector of coefficients of the canonical form
     */
    PolynomialFF construct_canonical(const uint zi, std::vector<PolynomialFF>& ai);
    /**
     *    Iterative construction of the canonical form
     *    @param zi the integer i to a z
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @param i The iteration step; stops at ai.size()
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(const uint zi, uint i, std::vector<PolynomialFF>& ai);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    /**
     *    Converts a vector of FFInts to a vector of mpz_class
     *    @param ai a vector of FFInts
     *    @return The vector ai converted to mpz_class objects
     */
    mpz_map convert_to_mpz(const PolynomialFF& poly) const;
    /**
     *    Convert a vector of RationalNumber objects to FFInts
     *    @param ri a vector of RationalNumber objects
     *    @param prime a prime defining the current finite field
     *    @return the vector ri converted to FFInt objects
     */
    ff_map convert_to_ffint(const rn_map& ri) const;
    /**
     * 
     */
    PolynomialFF solve_gauss();
    int deg = -1;
    uint n; /**< The number of parameters */
    bool use_chinese_remainder = false;
    bool check = false;
    Polynomial result;
    std::vector<std::vector<FFInt>> coef_mat {};
    std::vector<std::vector<uint>> rec_degs {};
    ff_map solved_degs {};
    mpz_class combined_prime; /**< The combination of the used prime numbers with the chinese remained theorem */
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    std::unordered_map<uint, std::vector<FFInt>> yis {};
    std::unordered_map<uint, std::vector<PolynomialFF>> ais {};
    std::unordered_map<uint, int> max_deg {};
  };
}
