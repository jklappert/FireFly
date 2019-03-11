// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "PolyReconst.hpp"
#include "RationalFunction.hpp"
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
    RatReconst(const RatReconst& other);
    RatReconst(RatReconst && other);
    RatReconst& operator=(const RatReconst& other);
    RatReconst& operator=(RatReconst && other);
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
  private:
    /**
     *  Starts the real interpolation managed by the class itself
     *  @param new_ti the currently used value of the homogenization variable
     *  @param num the black box probe
     *  @param feed_zi_ord the corresponding zi_order
     */
    void interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord);
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @param num f(y_i)
     *    @return a(i)
     */
    FFInt comp_ai(int i, int ip, const FFInt& num);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param rf the rational function
     *    @return A normalized version of ratFun
     */
    RationalFunction normalize(RationalFunction& rf);
    /**
     *    Constructs the canonical form of the rational function recursivly
     *    @return the rational function in its canonical form
     */
    std::pair<ff_map, ff_map>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param i an integer telling the current degree of the rational function
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical(uint32_t i);
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @returns f(y_i)
     */
    FFInt comp_fyi(uint32_t i, uint32_t ip, const FFInt& y);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
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
    std::vector<FFInt> ai {};
    std::unordered_map<uint32_t, PolyReconst> coef_n {};
    std::unordered_map<uint32_t, PolyReconst> coef_d {};
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_num {};// a vector entry should be just a pointer to save memory
    std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>> non_solved_degs_den {};
    polff_vec_map sub_num {};
    polff_vec_map sub_den {};
    ff_map_map saved_num_num {};
    ff_map_map saved_num_den {};
    int max_deg_num = -1;
    int max_deg_den = -1;
    int curr_deg_num = -1;
    int curr_deg_den = -1;
    std::vector<uint32_t> curr_zi_order_num {};
    std::vector<uint32_t> curr_zi_order_den {};
    uint32_t tmp_solved_coefs_num = 0;
    uint32_t tmp_solved_coefs_den = 0;
    void remove_ni(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
    void remove_di(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn);
    RationalFunction result;
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
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
     *  Solves a generalzied transposed Vandermonde system
     *  @param degs the monomial degrees of the system
     *  @param nums the probes of the polynomial
     *  @return a PolynomialFF object corresponding to the solution of the system
     */
    PolynomialFF solve_transposed_vandermonde(std::vector<std::vector<uint32_t>>& degs,
                                              const std::vector<std::pair<FFInt, uint32_t>>& nums);
    /**
     *  Solves a transposed Vandermonde system and calculates the terms originating from the shift after the first prime for the numerator
     *  @param key the current degree which has to be solved
     */
    void set_new_curr_deg_num_singular(uint32_t key);
    /**
     *  Solves a transposed Vandermonde system and calculates the terms originating from the shift after the first prime for the denominator
     *  @param key the current degree which has to be solved
     */
    void set_new_curr_deg_den_singular(uint32_t key);
    /**
     *  Saves the state of the current object and writes it to the specified file
     */
    void save_state();
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
    std::vector<bool> parsed_variables {std::vector<bool>(18, false)};
    int curr_parsed_variable = -1;
    uint32_t sub_count_num = 0;
    uint32_t sub_count_den = 0;
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
    bool scan = false;
    std::vector<uint32_t> all_shift_max_degs {};
    static std::vector<uint32_t> curr_shift;
    bool shift_works = false;
    bool normalize_to_den = true;
    int start_deg_num = 0;
    int start_deg_den = 1;
    uint32_t shifted_max_num_eqn = 0;
    std::unordered_set<uint32_t> shifted_degs_num {};
    std::unordered_set<uint32_t> shifted_degs_den {};
    std::unordered_set<uint32_t> zero_degs_num {};
    std::unordered_set<uint32_t> zero_degs_den {};
    bool normalizer_den_num = false;
    enum save_variables {COMBINED_PRIME, IS_DONE, MAX_DEG_NUM, MAX_DEG_DEN, NEED_PRIME_SHIFT,
                         NORMALIZER_DEG, NORMALIZE_TO_DEN, NORMALIZER_DEN_NUM, SHIFTED_MAX_NUM_EQN, SHIFT,
                         SHIFTED_DEGS_NUM, SHIFTED_DEGS_DEN, ZERO_DEGS_NUM, ZERO_DEGS_DEN, G_NI, G_DI,
                         COMBINED_NI, COMBINED_DI
                        };
  };
}
