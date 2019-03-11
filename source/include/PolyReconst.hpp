// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

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
  private:
    /**
     *  Starts the real interpolation managed by the class itself
     *  @param num the black box probe
     *  @param zi_ord the corresponding zi_order
     */
    void interpolate(const FFInt& num, const std::vector<uint32_t>& zi_ord);
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively
     *    @param tmp_zi the integer i to a zi
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @param num f(y_i)
     *    @param ai The vector of previously computed ai
     *    @return a(i)
     */
    PolynomialFF comp_ai(const uint32_t tmp_zi, int i, int ip, const PolynomialFF& num, std::vector<PolynomialFF>& ai);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param tmp_zi the integer i to a zi
     *    @param ai The computed ai
     *    @return The vector of coefficients of the canonical form
     */
    ff_map construct_canonical(const uint32_t tmp_zi, std::vector<PolynomialFF>& ai);
    /**
     *    Iterative construction of the canonical form
     *    @param tmp_zi the integer i to a z
     *    @param i The iteration step; stops at ai.size()
     *    @param ai The computed ai
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(const uint32_t tmp_zi, uint32_t i, std::vector<PolynomialFF>& ai);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param num a black box probe at a given parameter point
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    /**
     *  @return a PolynomialFF object of the solved transposed Vandermonde system
     */
    PolynomialFF solve_transposed_vandermonde();
    std::list<std::tuple<FFInt, std::vector<uint32_t>>> queue;
    int deg = -1;
    bool with_rat_reconst = false;
    Polynomial result;
    PolynomialFF result_ff;
    std::vector<std::vector<uint32_t>> rec_degs {};
    ff_map solved_degs {};
    std::vector<FFInt> nums {};
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    polff_vec_map ais {};
    std::unordered_map<uint32_t, int> max_deg {};
    static std::mutex mutex_statics;
    static ff_pair_map rand_zi;
  };
}
