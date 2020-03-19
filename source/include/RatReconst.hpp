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

#include "FFThieleInterpolator.hpp"
#include "PolyReconst.hpp"
#include "RationalFunction.hpp"
#include "RationalFunctionFF.hpp"

#ifdef FLINT
#include <flint/nmod_poly.h>
#endif

#include <map>

namespace firefly {
  /**
   * @class RatReconst
   * @brief A class to reconstruct a rational function from its values
   */
  class RatReconst : public BaseReconst {
  public:
    RatReconst() {}
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(uint32_t n_);
    /**
     *  Resets all static variables of the object
     *  @param change_prime when set to true, the prime of the FFInts is set to the first one of the prime array
     */
    static void reset(bool change_prime = true);
    /**
     *  Feeds a black-box probe which will be processed by the class.
     *  @param new_ti the value of t for the current feed
     *  @param num the black-box probe
     *  @param fed_zi_ord the corresponding zi_order to this feed
     *  @param fed_prime the corresponding prime number to this feed
     *  @return true if no interpolation is running, false otherwise, second true if one should write saved probes to file
     */
    std::pair<bool, bool> feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& fed_zi_ord, const uint32_t fed_prime);
    /**
     *  Feeds a black-box probe which will be processed by the class.
     *  @param new_ti the values of t for the current feed
     *  @param num the black-box probes
     *  @param fed_zi_ord the corresponding zi_orders to these feeds
     *  @param fed_prime the corresponding prime number to these feeds
     *  @return true if no interpolation is running, false otherwise, second true if one should write saved probes to file
     */
    std::pair<bool, bool> feed(const std::vector<FFInt>& new_ti, const std::vector<FFInt>& num, const std::vector<std::vector<uint32_t>>& fed_zi_ord, const uint32_t fed_prime);
    /**
     *  @return the result of the reconstruction as a RationalFunction object
     */
    RationalFunction get_result();
    /**
     *  @return the result of the rational function over the current field
     */
    RationalFunctionFF get_result_ff();
#ifdef FLINT
    /**
     *  @return the found factors over the current field, first entry is numerator, second is denominator
     */
    std::pair<std::vector<std::string>, std::vector<std::string>> get_factors_ff();
    /**
     *  @param factor_pos includes the positions of the factors that should be brought in canconical form, first entry is numerator, second is denominator
     *  @return Canonical form of real factors as a maps, where the first entry is the degree and the second is the coefficient for all non-zero entries
     */
    std::pair<std::unordered_map<uint32_t, uint64_t>, std::unordered_map<uint32_t, uint64_t>> get_canonical_factors(const std::pair<std::unordered_set<uint32_t>, std::unordered_set<uint32_t>>& factors_pos);
#endif
    /**
     *  Starts an interpolation job
     *  @return a tuple consisting of:
     *  bool true if it has done anything and false otherwise
     *  bool done
     *  uint32_t the prime counter
     */
    std::tuple<bool, bool, uint32_t> interpolate();
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
     *  @return the random number of zi at a given order
     */
    FFInt get_rand_zi(uint32_t zi, uint32_t order) const;
    /**
     *  @param order a vector of all zi orders
     *  @param generate generates the random values if they are not already generated
     *  @return a vector of all random numbers of the given orders
     */
    std::vector<FFInt> get_rand_zi_vec(const std::vector<uint32_t>& order, bool generate = false);
    /**
     *  @param zi the zi of which one wants to get the shift
     *  @return the corresponding shift as an FFInt
     */
    FFInt get_zi_shift(uint32_t zi) const;
    /**
     *  @return a vector which holds all shifts for each variable ordered like {z1_s, z2_s,...}, where zi_s is the shift of zi
     */
    std::vector<FFInt> get_zi_shift_vec() const;
    /**
     *  @param prime_counter indicates the prime number for which one checks if a shift is needed
     *  @return true if the reconstruction object still needs a shift
     */
    bool need_shift(uint32_t prime_counter);
    /**
     *  Sets a tag to save the state of the object after each prime
     *  @param tag_ the tag of the object
     */
    void set_tag(const std::string& tag_);
    /**
     *  Sets a user readable tag to this object
     *  @param tag_name_ the tag name of the object
     */
    void set_tag_name(const std::string& tag_name_);
    /**
     *  Get the tag of this object
     *  @return the tag of this object
     */
    std::string get_tag();
    /**
     *  Get the tag name of this object
     *  @return the tag name of this object
     */
    std::string get_tag_name();
    /**
     *  If called, the reconstruction will start from the state saved in a file
     *  @param file_name the absolute path to the saved state file
     *  @return a pair of a bool which indicates if the current objects needs a shift and an uint32_t for the prime number
     */
    std::pair<bool, uint32_t> start_from_saved_file(const std::string & file_name);
    /**
     *  Parses all saved probes for the current prime field
     *  @param file_name the absolute path to the probes save file
     *  @return a map of all t values feeded for a specific zi_order
     */
    std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> read_in_probes(const std::string& file_name);
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
     *  Sets the shift to the given tuple
     *  @param shift_ the given shift
     */
    void set_shift(const std::vector<FFInt>& shift_);
    /**
     *  @return is_interpolating
     */
    bool get_is_interpolating() const;
    /**
     *  With this option a full interpolation is performed over all prime fields to be sensitive of unlicky primes
     */
    void set_safe_interpolation();
    /**
     *  Returns the required feeds for new primes
     *  @return the vector of required feeds where the first entry is the multiplicity and the second entry is the number of different t, e.g, (3,4) means this object needs for 3 different zi order 4 black-box probes with different t
     */
    std::vector<std::pair<uint32_t, uint32_t>> get_needed_feed_vec();
    /**
     *  Returns the maximal degree of numerator and denominator
     *  @return The maximal degree of numerator (first) and denominator (second) as a pair
     */
    std::pair<uint32_t, uint32_t> get_max_deg();
    /**
     *  @return the anchor points
     */
    std::vector<FFInt> get_anchor_points();
    /**
     *  Sets the anchor points to the given tuple
     *  @param anchor_points the anchor points
     */
    void set_anchor_points(const std::vector<FFInt>& anchor_points);
    /**
     *  Writes current feed to file to reuse this information
     */
    void write_food_to_file();
    /**
     *  Calculates factors over the current prime field after the interpolation
     *  @param var_ sets the variable that should be replaced into the factors
     *  @param degs sets the occurring degs to avoid Thiele
     */
    void calc_factors(const std::string& var_ = "x", const std::pair<std::list<uint32_t>, std::list<uint32_t>>& degs = std::pair<std::list<uint32_t>, std::list<uint32_t>>());
    /**
     *  Sets the internal prime number to the maximum. This is just a hacky function for factor scans.
     */
    void set_prime_to_max();
    /**
     *  Returns a pair of a vector of pairs of a vector of required zi_orders and multiplicities and the maximum system size
     *  @return a pair of a vector of pairs of a vector of required zi_orders and multiplicities and the maximum system size
     */
    std::pair<std::vector<std::pair<std::vector<uint32_t>, uint32_t>>, uint32_t> get_zi_orders() const;
  private:
    /**
     *  Starts the real interpolation managed by the class itself
     *  @param new_ti the currently used value of the homogenization variable
     *  @param num the black box probe
     *  @param fed_zi_ord the corresponding zi_order
     */
    void interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& fed_zi_ord);
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
     *  @param rec the PolyReconst object which should be fed
     *  @param saved_num the container of saved probes for the polynomial
     *  @param is_num a bool which indicates if the polynomial is in the numerator (true) or denominator (false)
     *  @param sparse indicates whether the polynomial interpolation shoud be done sparsely
     *  @return true if the polymoial interpolation is done
     */
    bool feed_poly(uint32_t curr_deg, PolyReconst& rec, ff_map_map& saved_num, bool is_num, bool sparse = false);
    /**
     *  Combines the results of two different prime fields
     *  @param numerator the interpolated numerator of the current prime field
     *  @param denominator the interpolated denominator of the current prime field
     */
    void combine_primes(ff_map& numerator, ff_map& denominator);
    /**
     *  Builds a univariate system of equations for a rational function in the first prime only during factor scan
     *  @param tmp_ti the currently used value of t
     *  @param tmp_num the currently probed value of the black box
     */
    bool build_factor_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num);
    /**
     *  Solves the system of equations for the univariate factor rational function
     *  @return the solution of the system of equations
     */
    std::pair<ff_map, ff_map> solve_factor_uni_gauss();
    /**
     *  Builds a univariate system of equations for a rational function in the first prime
     *  @param tmp_ti the currently used value of t
     *  @param tmp_num the currently probed value of the black box
     *  @param yis the currently used tuple of yis
     */
    void build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis);
    /**
     *  Builds a univariate system of equations for a rational function for additional primes
     *  @param tmp_ti the currently used value of t
     *  @param tmp_num the currently probed value of the black box
     *  @param yis the currently used tuple of yis
     */
    void build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis);
    /**
     *  Adds a coefficient for the numerator which appears to be correctly reconstructed
     *  @param deg_vec the degree of the coefficient
     *  @param rn the coefficient
     */
    void remove_ni(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
    /**
     *  Adds a coefficient for the denominator which appears to be correctly reconstructed
     *  @param deg_vec the degree of the coefficient
     *  @param rn the coefficient
     */
    void remove_di(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
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
    /**
     *  Sets all required variables of this class in the "singular system" mode
     */
    void set_singular_system_vars();
    /**
     *  Checks if the reconstruction is done
     *  @param num a numerical value of the black box
     *  @param ti a numerical value for the probing variable
     *  @return true if the reconstruction is done
     */
    bool check_if_done(const FFInt& num, const FFInt& ti);
    /**
     *  Calculates the polynomials emerging from a parameter shift effecting lower degrees
     *  @param poly the seed polynomial
     *  @param deg the degree of the polynomial
     *  @return A map with a degree as the key and the polynomial emerging from the shift as its value
     */
    polff_map calculate_shift_polynomials(const PolynomialFF& poly, uint32_t deg);
    std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> queue; /**< A queue for all feeds */
    std::vector<std::vector<FFInt>> coef_mat {}; /**< A matrix that holds the system of equations for the homogenized system */
    std::unordered_map<uint32_t, std::vector<FFInt>> coef_mat_num {}; /**< A matrix that holds systems of equations for the Vandermonde systems of the numerator */
    std::unordered_map<uint32_t, std::vector<FFInt>> coef_mat_den {}; /**< A matrix that holds systems of equations for the Vandermonde systems of the denominator */
    PolynomialFF solved_num; /**< Stores already solved monomials of the numerator */
    PolynomialFF solved_den; /**< Stores already solved monomials of the denominator */
    ff_queue_map saved_ti {}; /**< Stores values of the homogenization variable t */
    std::list<std::pair<uint32_t, PolyReconst>> coef_n {}; /**< Stores PolyReconst objects for the numerator */
    std::list<std::pair<uint32_t, PolyReconst>> coef_d {}; /**< Stores PolyReconst objects for the denominator */
    // a vector entry should be just a pointer to save memory
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_num {}; /**< Stores unsolved degrees of the numerator */
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_den {}; /**< Stores unsolved degrees of the denominator */
    std::unordered_map<uint32_t, FFInt> num_sub_num {}; /**< Stores the numerical subtraction of the shift for each degree of the numerator */
    std::unordered_map<uint32_t, FFInt> num_sub_den {}; /**< Stores the numerical subtraction of the shift for each degree of the denominator */
    polff_map sub_num {}; /**< Stores the algebraic subtraction of the shift for each degree of the numerator */
    polff_map sub_den {}; /**< Stores the algebraic subtraction of the shift for each degree of the denominator */
    ff_map_map saved_num_num {}; /**< Stores interpolation points for each degree of the numerator */
    ff_map_map saved_num_den {}; /**< Stores interpolation points for each degree of the denominator */
    RationalFunction result; /**< Stores the result of the interpolation with rational numbers as coefficients */
    rn_map g_ni {}; /**< rational coefficient guesses for the numerator*/
    rn_map g_di {}; /**< rational coefficient guesses for the denominator*/
    mpz_map combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    mpz_map combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
    mpz_map combined_primes_ni {}; // used for safe mode
    mpz_map combined_primes_di {}; // used for safe mode
    polff_map solved_degs_num {}; /**< Stores solved polynomials of the numerator of the current field */
    polff_map solved_degs_den {}; /**< Stores solved polynomials of the denominator of the current field */
    ff_map tmp_res_num {}; /**< Stores the temporary numerator result interpolated over the current field */
    ff_map tmp_res_den {}; /**< Stores the temporary denominator result interpolated over the current field */
    std::vector<uint32_t> normalizer_deg {}; /**< Stores the degree which is used for normalization */
    std::string tag = ""; /**< The tag of this interpolation class for state saving */
    std::string tag_name = ""; /**< The tag name of this interpolation class for state saving */
    std::string var = ""; /**< The variable that is replaced when calculating factors */
    static std::vector<FFInt> shift; /**< The static shift used to obtain a unique normalization */
    static ff_pair_map rand_zi; /**< A static map storing potencies of the anchor points */
    static std::unordered_set<uint32_t> singular_system_set; /**< A static set which indicates if a shift is needed for a given prime field */
    static std::mutex mutex_statics;
    std::vector<bool> parsed_variables {std::vector<bool>(22, false)};  /**< A vector in which each entry indicates one variable which has to be parsed */
    std::vector<uint32_t> all_shift_max_degs {}; /**< Stores the maximal degree of numerator and denominator when shifting all variables */
    static std::vector<uint32_t> curr_shift; /**< Stores the current tested shift during a scan */
    std::pair<uint32_t, uint32_t> max_num_coef_num = std::make_pair(0, 0); // deg and number of terms
    std::pair<uint32_t, uint32_t> max_num_coef_den = std::make_pair(0, 0); // deg and number of terms
    std::unordered_set<uint32_t> dense_solve_degs_num {}; /**< Stores all degrees of the numerator which should be solved densely */
    std::unordered_set<uint32_t> dense_solve_degs_den {}; /**< Stores all degrees of the denominator which should be solved densely */
    std::unordered_set<uint32_t> zero_degs_num {}; /**< Stores all degrees of the numerator which have a zero coefficient */
    std::unordered_set<uint32_t> zero_degs_den {}; /**< Stores all degrees of the denominator which have a zero coefficient */
    std::vector<std::pair<uint32_t, uint32_t>> needed_feed_vec {}; /**< Stores all needed feeds for the next prime field */
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); /**< Timestamp that tracks the time until probes should be written to disk */
    ThieleInterpolator t_interpolator; /**< An object for Thiele interpolations */
    std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> saved_food; /**< Data structre used to write already used probes to a file from which one can resume if crashes occur. First FFInt is t second is num */
    std::map<std::vector<uint32_t>, std::vector<std::pair<uint64_t, uint64_t>>> parsed_probes {};
    std::pair<std::list<uint32_t>, std::list<uint32_t>> factor_degs {}; /**< Stores the degrees after one run with Thiele to build a system of equations instead */
#ifdef FLINT
    std::pair<std::vector<std::string>, std::vector<std::string>> factors; /**< Stores the found factors over the current field */
    std::pair<std::vector<std::pair<std::string, uint32_t>>, std::vector<std::pair<std::string, uint32_t>>> factors_flint; /**< Stores the found factors over the current field, second entry is exponent of factor */
#endif
    int max_deg_num = -1; /**< The maximal degree of the numerator */
    int max_deg_den = -1; /**< The maximal degree of the denominator */
    uint32_t vandermonde_size = 1; /**< Multiplicity of Vandermonde systems requried by PolyReconst */
    int curr_parsed_variable = -1;  /**< The current variable for the parser */
    uint32_t tmp_solved_coefs_num = 0; /**< A temporary variable that saves how many coefficients of the numerator have already been interpolated */
    uint32_t tmp_solved_coefs_den = 0; /**< A temporary variable that saves how many coefficients of the denominator have already been interpolated */
    uint32_t shifted_max_num_eqn = 0; /**< Stores the maximal number of equations when using a shift */
    uint32_t interpolations = 1;  /**< Indication how many interpolations should be made until one uses Vandermonde systems */
    size_t zero_counter = 0; /**< A counter for feeded zeros */
    bool is_singular_system = false; /**< Indicates whether one needs the shift to avoid singular systems */
    bool is_calc_factors = false; /**< Indicates whether factors should be calculted after the interpolation */
    bool from_save_state = false; /**< Indicates wether one resumes from a saved state */
    bool is_writing_probes = false; /**< Indicates wether this object is currently writing to a file */
    bool start_interpolation = true; /**< Indicats wether one can start the interpolation */
    bool first_run = true; /**< Indicates whether the current feed is the first interpolation run in the current prime field */
    bool scan = false; /**< Indicates whether this object scans for a sparse shift */
    bool shift_works = false; /**< Indicates whether the current shift works for normalization */
    bool normalize_to_den = true; /**< Indicates whether we normalize to the numerator or denominator */
    bool div_by_zero = false; /**< If the coefficient of the normalizer degree is zero, this variable forces the object to omit the prime field */
    bool first_feed = true; /**< Indicates whether this is the first feed in the current field */
    bool check_interpolation = false; /**< Indicates whether one has to check the interpolation point */
    bool fed_zero = false; /**< Indicates that the current feed was a zero */
    bool restart_sparse_interpolation_num = false; /**< Indicates whether one should proceed with a sparse interpolation instead of a dense one */
    bool restart_sparse_interpolation_den = false; /**< Indicates whether one should proceed with a sparse interpolation instead of a dense one */
    bool normalizer_den_num = false; /**< If true the real normalization degree is the denominator else the numerator */
    bool skip_thiele = false; /**< If true Thiele is skipped and replaced by system of equations during factor scan */
    enum save_variables {COMBINED_PRIME, TAG_NAME, IS_DONE, MAX_DEG_NUM, MAX_DEG_DEN, NEED_PRIME_SHIFT,
                         NORMALIZER_DEG, NORMALIZE_TO_DEN, NORMALIZER_DEN_NUM, SHIFTED_MAX_NUM_EQN, SHIFT,
                         SUB_NUM, SUB_DEN, ZERO_DEGS_NUM, ZERO_DEGS_DEN, G_NI, G_DI,
                         COMBINED_NI, COMBINED_DI, COMBINED_PRIMES_NI, COMBINED_PRIMES_DI, INTERPOLATIONS
                        };
  };
}
