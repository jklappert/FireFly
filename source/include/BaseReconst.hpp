#pragma once
#include <vector>
#include <mutex>
#include <unordered_map>
#include "RationalNumber.hpp"
#include "PolynomialFF.hpp"
#include "UintHasher.hpp"
#include "FFInt.hpp"

namespace firefly {
  typedef std::unordered_map<std::vector<uint>, mpz_class, UintHasher> mpz_map;
  typedef std::unordered_map<uint, std::unordered_map<std::vector<uint>, mpz_class, UintHasher>> mpz_map_map;
  typedef std::unordered_map<std::pair<uint, uint>, FFInt, UintPairHasher> ff_pair_map;
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;
  typedef std::unordered_map<std::vector<uint>, std::unordered_map<std::vector<uint>, FFInt, UintHasher>, UintHasher> ff_map_map;
  typedef std::unordered_map<std::vector<uint>, std::vector<std::pair<FFInt, FFInt>>, UintHasher> ff_vec_map;
  typedef std::unordered_map<std::vector<uint>, uint, UintHasher> uint_map;

  class BaseReconst {
  public:
    BaseReconst();
    BaseReconst(const BaseReconst& other);
    BaseReconst(BaseReconst && other);
    BaseReconst& operator=(const BaseReconst& other);
    BaseReconst& operator=(BaseReconst && other);
    FFInt get_rand();
    uint get_num_eqn();
    //void generate_anchor_points(uint max_order = 1);
    bool is_done();
    bool is_new_prime();
    uint get_prime();
    std::vector<uint> get_zi_order();
    uint get_zi();
  protected:
    bool use_chinese_remainder = false;
    bool check = false;
    bool done = false;
    bool new_prime = false;
    bool is_interpolating = false;
    std::vector<uint> curr_zi_order {};
    uint prime_number = 0;
    uint num_eqn = 0;
    uint n = 0; /**< The number of parameters */
    uint type;
    uint zi = 1;
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
  };
}
