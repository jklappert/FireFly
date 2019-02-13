#pragma once
#include <vector>
#include <mutex>
#include <unordered_map>
#include "RationalNumber.hpp"
#include "PolynomialFF.hpp"
#include "UintHasher.hpp"
#include "FFInt.hpp"
#include <stdint.h>

namespace firefly {
  typedef std::unordered_map<std::vector<uint32_t>, mpz_class, UintHasher> mpz_map;
  typedef std::unordered_map<uint32_t, std::unordered_map<std::vector<uint32_t>, mpz_class, UintHasher>> mpz_map_map;
  typedef std::unordered_map<std::pair<uint32_t, uint32_t>, FFInt, UintPairHasher> ff_pair_map;
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;
  typedef std::unordered_map<std::vector<uint32_t>, std::unordered_map<std::vector<uint32_t>, std::pair<FFInt, uint32_t>, UintHasher>, UintHasher> ff_map_map;
  typedef std::unordered_map<std::vector<uint32_t>, std::vector<std::pair<FFInt, FFInt>>, UintHasher> ff_vec_map;
  typedef std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> uint32_t_map;
  typedef std::unordered_map<uint32_t, std::vector<PolynomialFF>> polff_vec_map;
  typedef std::unordered_map<uint32_t, PolynomialFF> polff_map;

  class BaseReconst {
  public:
    BaseReconst();
    BaseReconst(const BaseReconst& other);
    BaseReconst(BaseReconst && other);
    BaseReconst& operator=(const BaseReconst& other);
    BaseReconst& operator=(BaseReconst && other);
    FFInt get_rand();
    uint32_t get_num_eqn();
    //void generate_anchor_points(uint32_t max_order = 1);
    bool is_done();
    bool is_new_prime();
    uint32_t get_prime();
    std::vector<uint32_t> get_zi_order();
    uint32_t get_zi();
  protected:
    bool use_chinese_remainder = false;
    bool check = false;
    bool done = false;
    bool new_prime = false;
    bool is_interpolating = false;
    std::vector<uint32_t> curr_zi_order {};
    uint32_t prime_number = 0;
    uint32_t num_eqn = 0;
    uint32_t n = 0; /**< The number of parameters */
    uint32_t type;
    uint32_t zi = 1;
    mpz_class combined_prime; /**< The combination of the used prime numbers with the chinese remained theorem */
    mutable std::mutex mutex_status;
    /**
     *    Converts the coefficients of a rational function from FFInts to mpz_class
     *    objects
     *    @param coefs a ff_map over a finite field
     *    @return the coefficients of the given rational function converted to
     *    mpz_class objects
     */
    mpz_map convert_to_mpz(const ff_map& coefs) const;
    /**
     *    Converts the elements of a vector of RationalNumber objects to FFInts
     *    @param ri the vector of RationalNumber objects
     *    @return elements of ri converted to FFInts
     */
    ff_map convert_to_ffint(const rn_map& ri) const;
    enum reconst_type {POLY, RAT};
  private:
    static std::mutex mutex_state;
    static uint64_t state;
    static uint64_t const multiplier;
    static uint64_t const increment;
    uint32_t rotr32(uint32_t x, uint32_t r);
    uint32_t pcg32();
    void pc32_init(uint64_t seed);
  };
}
