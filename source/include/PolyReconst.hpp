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

#pragma once

#include <cstdint>
#include <list>
#include "Polynomial.hpp"
#include "BaseReconst.hpp"

namespace firefly {
  /**
   * @class PolyReconst
   * @brief A class to reconstruct a polynomial from its values
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
     *    Feeds a new numerical value to the reconstruction algorithm
     *    @param yis the corresponding yi values to the feed
     *    @param num the numerical value of the black box
     */
    void feed(const std::vector<FFInt>& yis, const FFInt& num);
    /**
     *    Feeds a new numerical value to the reconstruction algorithm
     *    @param num the numerical value of the black box
     *    @param feed_zi_order the corresponding zi_order to the probe of the black box
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
    FFInt get_rand_zi(uint32_t zi, uint32_t order);
    /**
     *  @param orders a vector of all zi orders
     *  @return a vector of all random numbers of the given orders
     */
    std::vector<FFInt> get_rand_zi_vec(const std::vector<uint32_t>& orders);
    /**
     *  @return true if the rand_zi container is empty
     */
    bool is_rand_zi_empty();
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
    * Sets the threshold for the Ben-Or and Tiwari univariate interpolations
    * @param threshold the value of the threshold
    */
    static void set_bt_threshold(size_t threshold);
    /**
    * Turns the Newton Interpolation on/off. Default is on.
    * @param use_newton_new if true Newton Interpolation is used, if false it is not used
    */
    static void set_newton(bool use_newton_new);
    /**
    * Turns the Ben-Or and Tiwari Interpolation on/off. Default is on.
    * @param use_bt_new if true Newton Interpolation is used, if false it is not used
    */
    static void set_bt(bool use_bt_new);

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
     *    @param ip Recursion order
     *    @param num f(y_i)
     *    @param ai The vector of previously computed ai
     *    @return a(i)
     */
    FFInt comp_ai(int i, int ip, const FFInt& num, std::vector<FFInt>& ai);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param ai The computed ai
     *    @return The vector of coefficintients of the canonical form
     */
    ff_map construct_canonical(const std::vector<FFInt>& ai) const;
    /**
     *    Iterative construction of the canonical form
     *    @param i The iteration step; stops at ai.size()
     *    @param ai The computed ai
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(uint32_t i, const std::vector<FFInt>& ai) const;
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
    * updates the minimal polynomial with the Berlekamp-Massey algorithm for the new numerical
    * @param key the key of the corresponding Polynomials/Variables one wants to update
    * @return true if the termination condition is fulfilled and false otherwise
    */
    bool berlekamp_massey_step(std::vector<uint32_t>& key);
    /**
    * calculates the roots of the minimal polynomial (lambda) and the corresponding exponents to the anchor points
    * @param key the key of the Polynomial lambda on wich to find the roots on
    * @param base the base of wich the roots are powers of
    * @return a vector of the roots. They are stored as a pair where the first entry is the root itself and the second entry is the power to wich the anchor point equals aforementioned root
    */
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponents(std::vector<uint32_t>& key, const FFInt& base);
    /**
    * calculates the roots of the minimal polynomial (lambda) and the corresponding exponents to the anchor points, with the help of the Poly-Class; is not used, only exists for debugging
    * @param key the key of the Polynomial lambda on wich to find the roots on
    * @param base the base of wich the roots are powers of
    * @return a vector of the roots. They are stored as a pair where the first entry is the root itself and the second entry is the power to wich the anchor point equals aforementioned root
    */
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponents_with_poly_class(std::vector<uint32_t>& key, const FFInt& base);
    /**
    * solves the transposed shifted Vandermonde System
    * @param vis the evaluation point
    * @param fis the numerical values
    * @return a map with a degree as key and FFInt as value of the solved transposed Vandermonde system
    */
    std::vector<FFInt> solve_transposed_vandermonde(std::vector<FFInt>& vis, std::vector<FFInt>& fis);
    /**
    * if a Ben-Or and Tiwari interpolation suceeds, remove the degree from the vandermonde system
    * @param deg_vec the degree vector of the monomial of which the coefficient polynomial was successfully interpolated
    * @param coeffs the coefficients of the interpolated polynomial
    */
    void check_for_tmp_solved_degs_for_bt(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& coeffs, std::vector<size_t>& exponents);
    std::list<std::tuple<FFInt, std::vector<uint32_t>>> queue;
    int deg = -1;
    bool with_rat_reconst = false;
    Polynomial result;
    PolynomialFF result_ff;
    std::vector<std::vector<uint32_t>> rec_degs {};
    ff_map solved_degs {};
    ff_map tmp_solved_degs {};
    std::vector<FFInt> nums {};
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> ais {};
    std::unordered_map<uint32_t, int> max_deg {};
    static std::mutex mutex_statics;
    static ff_pair_map rand_zi;
    std::vector<uint32_t> zero_element {};
    bool combine_res = false;
    ff_map construct_tmp_canonical(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) const;
    void check_for_tmp_solved_degs_for_newton(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai);
    static size_t bt_threshold; // the amount of times the generator polynomial has to remain unchanged when the berlekamp_massey_step function gets called to terminate the Berlekamp/Massey algorithm, default is 1
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher> bt_terminator; // the amount of times the generator polynomial has not changed though the berlekamp_massey_step function got called
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher>  b; // the generator polynomial before it changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher>  l; // the amount of times the berlekamp_massey_step method got called ince the generator polynomial changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, FFInt, UintHasher>  delta; // the discrepancy when the generator polynomial changed its degree the last time
    std::unordered_map<std::vector<uint32_t>, size_t, UintHasher>  bm_iteration; // the amount of times the berlekamp_massey_step method has been called for a specific key
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> nums_for_bt; // stores the numerical values for the univariate interpolations used in Ben-Or and Tiwari
    std::unordered_map<std::vector<uint32_t>, std::vector<FFInt>, UintHasher> lambda; // The current generator polynomials in the Berlekamp-Massey algorithm
    static bool use_bt; // determines, if Ben-Or and Tiwari is used for univariate interpolation, default is true
    static bool use_newton; // determines, if Newton is used for univariate interpolation, default is true, if rational reconstruction is used this has always to be true
  };
}
