#pragma once

#include <cstdint>
#include <list>
#include "Polynomial.hpp"
#include "BaseReconst.hpp"

namespace firefly {
  class PolyReconst : public BaseReconst {
  public:
    /**
     *    A constructor
     *    @param n_ The number of parameters as an integer
     */
    PolyReconst(uint n_, const int deg_inp = -1, const bool with_rat_reconst_inp = false);
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
    void feed(const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime);
    void interpolate();
    FFInt get_rand_zi(uint zi, uint order);
    bool is_rand_zi_empty();
    Polynomial get_result();
    PolynomialFF get_result_ff();
    void generate_anchor_points();
    void set_anchor_points(const std::vector<FFInt> &anchor_points, bool force = false);
  private:
    void interpolate(const FFInt& num, const std::vector<uint>& zi_ord);
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param tmp_zi the integer i to a zi
     *    @param ai The vector of previously computed ai
     *    @param num f(y_i)
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @return a(i)
     */
    PolynomialFF comp_ai(const uint tmp_zi, int i, int ip, const PolynomialFF& num, std::vector<PolynomialFF>& ai);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param tmp_zi the integer i to a zi
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @return The vector of coefficients of the canonical form
     */
    ff_map construct_canonical(const uint tmp_zi, std::vector<PolynomialFF>& ai);
    /**
     *    Iterative construction of the canonical form
     *    @param tmp_zi the integer i to a z
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @param i The iteration step; stops at ai.size()
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(const uint tmp_zi, uint i, std::vector<PolynomialFF>& ai);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    /**
     * 
     */
    PolynomialFF solve_transposed_vandermonde();
    std::list<std::tuple<FFInt, std::vector<uint>>> queue;
    int deg = -1;
    bool with_rat_reconst = false;
    Polynomial result;
    PolynomialFF result_ff;
    std::vector<std::vector<uint>> rec_degs {};
    ff_map solved_degs {};
    std::vector<FFInt> nums {};
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    std::unordered_map<uint, std::vector<PolynomialFF>> ais {};
    std::unordered_map<uint, int> max_deg {};
    static std::mutex mutex_statics;
    static ff_pair_map rand_zi;
  };
}
