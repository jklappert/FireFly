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

#include "PolyReconst.hpp"
#include "RationalFunction.hpp"
#include "FFThieleInterpolator.hpp"
#include <unordered_set>
#include <queue>

namespace firefly {
  /**
   * @class RatReconst
   * @brief A class to reconstruct a rational function from its values
   */
  class RatReconst : public BaseReconst {
  public:
    RatReconst() {};
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(uint32_t n_);
    /**
     *  Resets all static variables of the object
     */
    static void reset();
    void feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t fed_prime);
    /**
     *  @return the result of the reconstruction as a RationalFunction object
     */
    RationalFunction get_result();
    /**
     *  Starts an interpolation job
     */
    bool interpolate();
    /**
     *  Disables the shift, thus setting it to a zero vector
     */
    void disable_shift();
    /**
     *  Generates anchor points
     */
    void generate_anchor_points();
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
    std::vector<FFInt> get_rand_zi_vec(const std::vector<uint32_t>& order);
    /**
     *  @param zi the zi of which one wants to get the shift
     *  @return the corresponding shift as an FFInt
     */
    FFInt get_zi_shift(uint32_t zi);
    /**
     *  @return a vector which holds all shifts for each variable ordered like {z1_s, z2_s,...}, where zi_s is the shift of zi
     */
    std::vector<FFInt> get_zi_shift_vec();
    /**
     *  @return true if the reconstruction object still needs a shift
     */
    bool need_shift();
    /**
     *  Sets a tag to save the state of the object after each prime
     *  @param tag_ the tag of the object
     */
    void set_tag(const std::string& tag_);
    /**
     *  If called, the reconstruction will start from the state saved in a file
     *  @param file_name the absolute path to the saved state file
     */
    void start_from_saved_file(std::string file_name);
    /**
     *  Enables the scan for a sparsest shift
     */
    void scan_for_sparsest_shift();
    /**
     *  Sets a shift of a given n-tuple of 0 and 1
     *  @param shifted_zis a vector consisting of 0 and 1 corresponding to the variables which should be shifted
     */
    void set_zi_shift(const std::vector<uint32_t>& shifted_zis);
    /**
     *  @return true if the currently used shift produces a constant in the numerator or denominator
     */
    bool is_shift_working();
    /**
     *  Sets the currently set shift to be the one used during the reconstruction
     */
    void accept_shift();
    /**
     *  @return is_interpolating
     */
    bool get_is_interpolating();
    /**
     *  With this option a full interpolation is performed over all prime fields to be sensitive of unlicky primes
     */
    void set_safe_interpolation();
  private:
    /**
     *  Starts the real interpolation managed by the class itself
     *  @param new_ti the currently used value of the homogenization variable
     *  @param num the black box probe
     *  @param feed_zi_ord the corresponding zi_order
     */
    void interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param rf the rational function
     *    @return A normalized version of ratFun
     */
    RationalFunction normalize(RationalFunction& rf);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param num the numerical value of the black bock
     *    @param ti the homogenization variable t
     *    @return true or false
     */
    bool test_guess(const FFInt& num, const FFInt& ti);
    /**
     *  @return true if all coefficients can be promoted to Q
     */
    bool rec_rat_coef();
    /**
     *  Solves a univariate system of equations (only used for first prime)
     *  @return the solved coefficients ordered as numerator and denominator
     */
    std::pair<ff_map, ff_map> solve_gauss();
    /**
     *  Solves a univariate system of equations (only used for additional primes)
     *  @return the solved coefficients ordered as numerator and denominator
     */
    std::pair<ff_map, ff_map> solve_homogenized_multi_gauss();
    /**
     *  Feeds a PolyReconst object with the results of the univariate system of equations
     *  @param curr_deg the degree of the polynomial
     *  @param max_deg the maximal degree of the numerator (denominator)
     *  @param coef the map of corresponding PolyReconst objects
     *  @param rec the PolyReconst object which should be fed
     *  @param saved_num the container of saved probes for the polynomial
     *  @param sub_save the container of subtraction terms due to the shift
     *  @param is_num a bool which indicates if the polynomial is in the numerator (true) or denominator (false)
     */
    std::tuple<int, uint32_t, std::vector<uint32_t>> feed_poly(int curr_deg,
                                                               uint32_t max_deg, std::unordered_map<uint32_t, PolyReconst>& coef,
                                                               PolyReconst& rec, ff_map_map& saved_num, polff_vec_map& sub_save, bool is_num);
    void combine_primes(ff_map& numerator, ff_map& denominator);
    /**
     *  Builds a univariate system of equations for a rational function in the first prime
     *  @param tmp_ti the currently used value of t
     *  @param tmp_num the currently probed value of the black box
     *  @param yis the currently used tuple of yis
     */
    void build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis);
    /**
     *  Builds a univariate system of equations for a rational function for additional primes
     *  @param tmp_ti the currently used value of t
     *  @param tmp_num the currently probed value of the black box
     *  @param yis the currently used tuple of yis
     */
    void build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis);
    bool first_run = true;
    std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> queue;
    std::vector<std::vector<FFInt>> coef_mat {};
    std::unordered_map<uint32_t, std::vector<std::pair<FFInt, uint32_t>>> coef_mat_num {};
    std::unordered_map<uint32_t, std::vector<std::pair<FFInt, uint32_t>>> coef_mat_den {};
    PolynomialFF solved_num;
    PolynomialFF solved_den;
    uint32_t curr_zi = 2;
    ff_vec_map saved_ti {};
    std::unordered_map<uint32_t, PolyReconst> coef_n {};
    std::unordered_map<uint32_t, PolyReconst> coef_d {};
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_num {};// a vector entry should be just a pointer to save memory
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_den {};
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_num_copy {};// a vector entry should be just a pointer to save memory
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_den_copy {};
    std::unordered_map<uint32_t, FFInt> num_sub_num {};
    std::unordered_map<uint32_t, FFInt> num_sub_den {};
    polff_vec_map sub_num {};
    polff_vec_map sub_den {};
    ff_map_map saved_num_num {};
    ff_map_map saved_num_den {};
    int max_deg_num = -1;
    int max_deg_den = -1;
    int curr_deg_num = -1;
    int curr_deg_den = -1;
    uint32_t tmp_sol_const_num = 0;
    uint32_t tmp_sol_const_den = 0;
    std::vector<uint32_t> curr_zi_order_num {};
    std::vector<uint32_t> curr_zi_order_den {};
    FFInt shifted_const = 0;
    bool remove_const = false;
    uint32_t tmp_solved_coefs_num = 0;
    uint32_t tmp_solved_coefs_den = 0;
    uint32_t sub_count_num = 0;
    uint32_t sub_count_den = 0;
    void remove_ni(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
    void remove_di(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
    RationalFunction result;
    rn_map g_ni {}; /**< rational coefficient guesses for the numerator*/
    rn_map g_di {}; /**< rational coefficient guesses for the denominator*/
    mpz_map combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    mpz_map combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
    static std::mutex mutex_statics;
    /**
     *  Adds non-solved monomials of the numerator to a data object
     *  @param deg the monomial degree
     */
    void add_non_solved_num(const std::vector<uint32_t>& deg);
    /**
     *  Adds non-solved monomials of the denominator to a data object
     *  @param deg the monomial degree
     */
    void add_non_solved_den(const std::vector<uint32_t>& deg);
    /**
     *  Checks if complete univariate mapped degrees are solved to possibly omit the shift
     *  @param uni_degs a vector of univariate mapped degress
     *  @param is_num indicates if this is the numerator (true) or denominator (false)
     */
    void check_for_solved_degs(const std::vector<uint32_t>& uni_degs, const bool is_num);
    /**
     *  Saves the state of the current object and writes it to the specified file
     */
    void save_state();
    /**
     *  Saves the state of the current object when it encountered a zero in a consecutive prime field
     */
    void save_zero_consecutive_prime();
    /**
     *  Saves the state of the current object when it encountered a zero in the first prime field
     */
    void save_zero_state();
    polff_map solved_degs_num {};
    polff_map solved_degs_den {};
    std::vector<uint32_t> normalizer_deg {};
    FFInt const_den;
    std::string tag = "";
    bool is_singular_system = false;
    static std::vector<FFInt> shift;
    static ff_pair_map rand_zi;
    static bool need_prime_shift;
    static bool set_singular_system;
    void set_singular_system_vars();
    std::vector<bool> parsed_variables {std::vector<bool>(19, false)};
    int curr_parsed_variable = -1;
    /**
     *  Parses a vector from a file with a given number of maximal entries
     *  @param line a string representing the line which should be parsed
     *  @param number_of_parameters a limiting number how many entries should be parsed
     *  @return the parsed vector
     */
    std::vector<uint32_t> parse_vector(std::string& line, int number_of_parameters = -1);
    /**
     *  Parses a rational number from a file
     *  @param line the string that should be parsed to a rational number
     */
    std::vector<mpz_class> parse_rational_number(std::string& line);
    /**
     *  Parses a prime number counter from a file
     *  @param file_name the file name
     */
    void parse_prime_number(std::string& file_name);
    /**
     *  Checks if the reconstruction is done
     *  @param num a numerical value of the black box
     *  @param ti a numerical value for the probing variable
     *  @return true if the reconstruction is done
     */
    bool check_if_done(const FFInt& num, const FFInt& ti);
    /**
     *  @returns the anchor points
     */
    std::vector<FFInt> get_anchor_points();
    bool scan = false;
    std::vector<uint32_t> all_shift_max_degs {};
    static std::vector<uint32_t> curr_shift;
    bool shift_works = false;
    bool normalize_to_den = true;
    int start_deg_num = 0;
    int start_deg_den = 1;
    uint32_t shifted_max_num_eqn = 0;
    bool div_by_zero = false;
    bool first_feed = true;
    size_t zero_counter = 0;
    bool check_interpolation = false;
    bool is_zero = false;
    bool fed_zero = false;
    std::pair<uint32_t, uint32_t> max_num_coef_num = std::make_pair(0, 0); // deg and number of terms
    std::pair<uint32_t, uint32_t> max_num_coef_den = std::make_pair(0, 0); // deg and number of terms
    std::unordered_set<uint32_t> dense_solve_degs_num {};
    std::unordered_set<uint32_t> dense_solve_degs_den {};
    std::unordered_set<uint32_t> shifted_degs_num {};
    std::unordered_set<uint32_t> shifted_degs_den {};
    std::unordered_set<uint32_t> zero_degs_num {};
    std::unordered_set<uint32_t> zero_degs_den {};
    /**
     *  Calculates the polynomials emerging from a parameter shift effecting lower degrees
     *  @param poly the seed polynomial
     *  @param deg the degree of the polynomial
     *  @return A map with a degree as the key and the polynomial emerging from the shift as its value
     */
    polff_map calculate_shift_polynomials(const PolynomialFF& poly, uint32_t deg);
    bool normalizer_den_num = false;
    ThieleInterpolator t_interpolator;
    uint32_t interpolations = 1;
    enum save_variables {COMBINED_PRIME, IS_DONE, MAX_DEG_NUM, MAX_DEG_DEN, NEED_PRIME_SHIFT,
                         NORMALIZER_DEG, NORMALIZE_TO_DEN, NORMALIZER_DEN_NUM, SHIFTED_MAX_NUM_EQN, SHIFT,
                         SHIFTED_DEGS_NUM, SHIFTED_DEGS_DEN, ZERO_DEGS_NUM, ZERO_DEGS_DEN, G_NI, G_DI,
                         COMBINED_NI, COMBINED_DI, INTERPOLATIONS
                        };
  };
}
