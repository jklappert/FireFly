//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#pragma once

#include "firefly/BaseReconst.hpp"
#include "firefly/Polynomial.hpp"

#include <list>

namespace firefly {
  /**
   *  @class PolyReconst
   *  @brief A class to reconstruct a polynomial from its values. Note that this class only works as shown in RatReconst, i.e. set anchor points on a dummy instance before using this class for interpolation
   */
  class PolyReconst : public BaseReconst {
  public:
    /**
     *    A constructor
     *    @param n_ The number of parameters as an integer
     *    @param deg_inp the expected maximal degree of the black box
     */
    PolyReconst(uint32_t n_, const int deg_inp = -1, const bool with_rat_reconst_inp = false);
    /**
     *    Default constructor. Should not be used explicitly.
     */
    PolyReconst();
    /**
     *    Feeds a new numerical value to the reconstruction algorithm. Only useful in combination with RatRconst
     *    @param num the numerical value of the black box
     */
    void feed(const FFInt& num);
    /**
     *    Feeds a new numerical value to the reconstruction algorithm
     *    @param num the numerical value of the black box
     *    @param feed_zi_ord the corresponding zi_order to the probe of the black box
     *    @param fed_prime the counter of the prime number corresponding to the feed
     */
    void feed(const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t fed_prime);
    /**
     *  Starts an interpolation job
     */
    void interpolate();
    /**
     *  @param zi the zi of which one wants to get the corresponding random number
     *  @param order the order of zi, i.e. zi^order
     *  @returns the random number of zi at a given order
     */
    FFInt get_rand_zi(uint32_t zi, uint32_t order) const;
    /**
     *  @param orders a vector of all zi orders
     *  @return a vector of all random numbers of the given orders
     */
    std::vector<FFInt> get_rand_zi_vec(const std::vector<uint32_t>& orders) const;
    /**
     *  @return true if the rand_zi container is empty
     */
    bool is_rand_zi_empty() const;
    /**
     *  @return a Polynomial object if the reconstruction finished succesfully
     */
    Polynomial get_result();
    /**
     *  @return a PolynomialFF object at any time during the reconstruction. Not needed for the user.
     */
    PolynomialFF get_result_ff();
    /**
     *  Generates new anchor points
     */
    void generate_anchor_points();
    /**
     *  Sets the anchor points to given values
     *  @param anchor_points the values to which the anchor points should be set
     */
    void set_anchor_points(const std::vector<FFInt>& anchor_points, bool force = false);
    /**
     *  Resets all statics
     */
    static void reset();
    /**
     *  Sets the threshold for the Ben-Or and Tiwari univariate interpolations
     *  @param threshold the value of the threshold
     */
    void set_bt_threshold(size_t threshold);
    /**
     *  Turns the Newton Interpolation on/off. Default is on.
     *  @param use_newton_new if true Newton Interpolation is used, if false it is not used
     */
    void set_newton(bool use_newton_new);
    /**
     *  Turns the Ben-Or and Tiwari Interpolation on/off. Default is on.
     *  @param use_bt_new if true Newton Interpolation is used, if false it is not used
     */
    void set_bt(bool use_bt_new);
    /**
     *  Returns how many additional equations are required to solve the Vandermonde system
     *  @return the required number of equations
     */
    uint32_t get_vandermonde_num_eqn() const;
    /**
     *  Sets individual degree bounds on each variable
     *  @param individual_degree_bounds_ the degree bounds 
     */
    void set_individual_degree_bounds(const std::vector<uint32_t>& individual_degree_bounds_);
  private:
    /**
     *  Starts the real interpolation managed by the class itself
     *  @param num the black box probe
     *  @param zi_ord the corresponding zi_order
     */
    void interpolate(const FFInt& num, const std::vector<uint32_t>& zi_ord);
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively
     *    @param i The order of a(i)
     *    @param num f(y_i)
     *    @param ai The vector of previously computed ai
     *    @return a(i)
     */
    FFInt comp_ai(int i, const FFInt& num, const std::vector<FFInt>& ai) const;
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param ai The computed ai
     *    @return The vector of coefficintients of the canonical form
     */
    ff_map construct_canonical(const std::vector<FFInt>& ai) const;
    /**
     *    Iterative construction of the canonical form
     *    @param ai The computed ai
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(const std::vector<FFInt>& ai) const;
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param num a black box probe at a given parameter point
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    /**
     *  @return a map with a degree as key and FFInt as value of the solved transposed Vandermonde system
     */
    ff_map build_and_solve_transposed_vandermonde();
    /**
     *  updates the minimal polynomial with the Berlekamp-Massey algorithm for the new numerical
     *  @param key the key of the corresponding Polynomials/Variables one wants to update
     *  @return true if the termination condition is fulfilled and false otherwise
     */
    bool berlekamp_massey_step(const std::vector<uint32_t>& key);
    /**
     *  calculates the roots of the minimal polynomial (lambda) and the corresponding exponents to the anchor points
     *  @param key the key of the Polynomial lambda on wich to find the roots on
     *  @param base the base of wich the roots are powers of
     *  @return a vector of the roots. They are stored as a pair where the first entry is the root itself and the second entry is the power to wich the anchor point equals aforementioned root
     */
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponents(const std::vector<uint32_t>& key, const FFInt& base);
    /**
     *  calculates the roots of the minimal polynomial (lambda) and the corresponding exponents to the anchor points, with the help of the Poly-Class; is not used, only exists for debugging
     *  @param key the key of the Polynomial lambda on wich to find the roots on
     *  @param base the base of wich the roots are powers of
     *  @return a vector of the roots. They are stored as a pair where the first entry is the root itself and the second entry is the power to wich the anchor point equals aforementioned root
     */
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponents_with_poly_class(const std::vector<uint32_t>& key, const FFInt& base);
    /**
     *  solves the transposed shifted Vandermonde System
     *  @param vis the evaluation point
     *  @param fis the numerical values
     *  @return a map with a degree as key and FFInt as value of the solved transposed Vandermonde system
     */
    std::vector<FFInt> solve_transposed_vandermonde(const std::vector<FFInt>& vis, const std::vector<FFInt>& fis) const;
    /**
     *  if a Ben-Or and Tiwari interpolation suceeds, remove the degree from the vandermonde system
     *  @param deg_vec the degree vector of the monomial of which the coefficient polynomial was successfully interpolated
     *  @param coeffs the coefficients of the interpolated polynomial
     */
    void check_for_tmp_solved_degs_for_bt(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& coeffs, const std::vector<size_t>& exponents);
    ff_map construct_tmp_canonical(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) const;
    void check_for_tmp_solved_degs_for_newton(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai);
    std::list<std::tuple<FFInt, std::vector<uint32_t>>> queue; /**< A queue which holds all unfeeded feeds */
    Polynomial result; /**< The result of the interpolation with rational number coefficients */
    PolynomialFF result_ff; /**< The result of the interpolation with coefficients in the current field */
    std::vector<std::vector<uint32_t>> rec_degs {};
    ff_map solved_degs {}; /**< A map holding already solved degrees used for permanent pruning*/
    ff_map tmp_solved_degs {}; /**< A map holding degrees used for temporary pruning */
    std::vector<FFInt> nums {}; /**< A vector holding evaluations of the polynomial */
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> ais {}; /**<< Coefficients used for Newton interpolation */
    std::unordered_map<uint32_t, int> max_deg {};
    static std::mutex mutex_anchor; /**< A mutex for writing to global_anchor_points */
    static std::vector<FFInt> global_anchor_points;
    std::vector<FFInt> private_anchor_points;
    std::vector<uint32_t> zero_element {};
    size_t bt_threshold = 1; // the amount of times the generator polynomial has to remain unchanged when the berlekamp_massey_step function gets called to terminate the Berlekamp/Massey algorithm, default is 1
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher> bt_terminator; // the amount of times the generator polynomial has not changed though the berlekamp_massey_step function got called
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher>  b; // the generator polynomial before it changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher>  l; // the amount of times the berlekamp_massey_step method got called ince the generator polynomial changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, FFInt, UintHasher>  delta; // the discrepancy when the generator polynomial changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher>  bm_iteration; // the amount of times the berlekamp_massey_step method has been called for a specific key
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> nums_for_bt; // stores the numerical values for the univariate interpolations used in Ben-Or and Tiwari
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> lambda; // The current generator polynomials in the Berlekamp-Massey algorithm
    int deg = -1; /**< The maximal degree of the to be interpolated polynomial */
    bool with_rat_reconst = false; /**< A variable indicating whether this object is called from a RatReconst object */
    bool combine_res = false; /**< A bool indicating whether the result should be combined using the Chinese Remainder Theorem */
    bool use_bt = true; // determines, if Ben-Or and Tiwari is used for univariate interpolation, default is true
    bool use_newton = true; // determines, if Newton is used for univariate interpolation, default is true, if rational reconstruction is used this has always to be true
    bool is_set_individual_degree_bounds = false; /**< If true individual degree bounds are set */
    std::vector<uint32_t> individual_degree_bounds {}; /**< Stores individual degree bounds */
  };
}
