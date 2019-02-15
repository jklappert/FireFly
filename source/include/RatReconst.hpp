#pragma once

#include "PolyReconst.hpp"
#include "RationalFunction.hpp"
#include <unordered_set>

namespace firefly {

  typedef std::unordered_map<std::vector<uint32_t>, std::unordered_map<uint32_t, std::unordered_map<uint32_t, FFInt>>, UintHasher> shift_map;

  class RatReconst : public BaseReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(uint32_t n_);
    RatReconst(const RatReconst& other);
    RatReconst(RatReconst && other);
    RatReconst& operator=(const RatReconst& other);
    RatReconst& operator=(RatReconst && other);
    void feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t& fed_prime);
    RationalFunction get_result();
    void interpolate();
    void disable_shift();
    void generate_anchor_points();
    FFInt get_rand_zi(uint32_t zi, uint32_t order);
    std::vector<FFInt> get_rand_zi_vec(std::vector<uint32_t> order);
    FFInt get_zi_shift(uint32_t zi);
    std::vector<FFInt> get_zi_shift_vec();
    bool need_shift();
    void set_tag(std::string tag_);
    void start_from_saved_file(std::string file_name);
  private:
    void interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord);
    FFInt comp_ai(int i, int ip, const FFInt& num);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param ratFun the rational function
     *    @param prime a prime number defining the finite field
     *    @return A normalized version of ratFun
     */
    RationalFunction normalize(RationalFunction& rf);
    /**
     *    Constructs the canonical form of the rational function recursivly
     *    @param ai a vector of the coefficients ai as FFints
     *    @param prime a prime number defining the current finite field
     *    @return the rational function in its canonical form
     */
    std::pair<ff_map, ff_map>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param ai a vector of the coefficients ai as FFInts
     *    @param i an integer telling the current degree of the rational function
     *    @param prime a prime number defining the current finite field
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical(uint32_t i);
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param ai a vector of coefficients ai
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @param prime a prime number defining the finite field
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
    bool rec_rat_coef();
    std::pair<ff_map, ff_map> solve_gauss();
    std::pair<ff_map, ff_map> solve_homogenized_multi_gauss();
    std::tuple<int, uint32_t, std::vector<uint32_t>> feed_poly(int curr_deg,
                                                       uint32_t max_deg, std::unordered_map<uint32_t, PolyReconst>& coef,
                                                       PolyReconst& rec, ff_map_map& saved_num, polff_vec_map& sub_save, bool is_num);
    void combine_primes(ff_map& numerator, ff_map& denominator);
    void build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis);
    void build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis);
    bool first_run = true;
    std::list<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> queue;
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
    void remove_ni(const std::vector<uint32_t>& deg_vec, RationalNumber& rn);
    void remove_di(const std::vector<uint32_t>& deg_vec, RationalNumber& rn);
    RationalFunction result;
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
    rn_map g_ni {}; /**< rational coefficient guesses for the numerator*/
    rn_map g_di {}; /**< rational coefficient guesses for the denominator*/
    mpz_map combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    mpz_map combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
    static std::mutex mutex_statics;
    void add_non_solved_num(const std::vector<uint32_t>& deg);
    void add_non_solved_den(const std::vector<uint32_t>& deg);
    void check_for_solved_degs(std::vector<uint32_t>& uni_degs, const bool is_num);
    PolynomialFF solve_transposed_vandermonde(std::vector<std::vector<uint32_t>>& degs,
                                              const std::vector<std::pair<FFInt, uint32_t>>& nums);
    std::vector<FFInt> solve_uni_transposed_vandermonde(const std::vector<FFInt>& nums);
    void set_new_curr_deg_num_singular(uint32_t key);
    void set_new_curr_deg_den_singular(uint32_t key);
    polff_map solved_degs_num {};
    polff_map solved_degs_den {};
    std::vector<uint32_t> min_deg_den_vec {};
    FFInt const_den = 0;
    std::string tag = "";
    bool is_singular_system = false;
    static std::vector<FFInt> shift;
    static ff_pair_map rand_zi;
    static bool need_prime_shift;
    static bool set_singular_system;
    void set_singular_system_vars();
    std::vector<bool> parsed_variables {std::vector<bool>(9, false)};
    int curr_parsed_variable = -1;//new
    uint32_t sub_count_num = 0;//new
    uint32_t sub_count_den = 0;//new
    std::vector<uint32_t> parse_vector(std::string& line, int number_of_parameters = -1);
    std::vector<mpz_class> parse_rational_number(std::string& line);
    void parse_prime_number(std::string& file_name);
    /*shift_map saved_shifts_num {}; //new
    shift_map saved_shifts_den {}; //new
    std::unordered_set<uint32_t> zero_degs_num {}; //new
    std::unordered_set<uint32_t> zero_degs_den {}; //new
    uint32_t sub_count_num = 0;//new
    uint32_t sub_count_den = 0;//new*/
    //FFInt get_particular_shift(const std::vector<uint32_t>& zi_order, int deg, bool is_num, uint32_t sub_count);
    void calculate_shift(const PolynomialFF& poly, const std::vector<uint32_t>& zi_order, int deg, bool is_num);
    enum save_variables {COMBINED_PRIME, MAX_DEG_NUM, MAX_DEG_DEN, NEED_PRIME_SHIFT,
    MIN_DEG_DEN_VEC, G_NI, G_DI, COMBINED_NI, COMBINED_DI};
  };
}
