#include "BaseReconst.hpp"
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include <random>

namespace firefly {
  uint64_t const BaseReconst::multiplier = 6364136223846793005u;
  uint64_t const BaseReconst::increment  = 1442695040888963407u;
  uint64_t BaseReconst::state = 0x4d595df4d0f33173;
  std::mutex BaseReconst::mutex_state;

  BaseReconst::BaseReconst() {}

  BaseReconst::BaseReconst(const BaseReconst& other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    done = other.done;
    new_prime = other.new_prime;
    check = other.check;
    use_chinese_remainder = other.use_chinese_remainder;
    curr_zi_order = other.curr_zi_order;
    prime_number = other.prime_number;
    num_eqn = other.num_eqn;
    n = other.n;
    type = other.type;
    zi = other.zi;
    combined_prime = other.combined_prime;
  }

  BaseReconst::BaseReconst(BaseReconst && other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    done = std::move(other.done);
    new_prime = std::move(other.new_prime);
    check = std::move(other.check);
    use_chinese_remainder = std::move(other.use_chinese_remainder);
    curr_zi_order = std::move(other.curr_zi_order);
    prime_number = std::move(other.prime_number);
    num_eqn = std::move(other.num_eqn);
    n = std::move(other.n);
    type = std::move(other.type);
    zi = std::move(other.zi);
    combined_prime = std::move(other.combined_prime);
  }

  BaseReconst& BaseReconst::operator=(const BaseReconst& other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      done = other.done;
      new_prime = other.new_prime;
      check = other.check;
      use_chinese_remainder = other.use_chinese_remainder;
      curr_zi_order = other.curr_zi_order;
      prime_number = other.prime_number;
      num_eqn = other.num_eqn;
      n = other.n;
      type = other.type;
      zi = other.zi;
      combined_prime = other.combined_prime;
    }

    return *this;
  }

  BaseReconst& BaseReconst::operator=(BaseReconst && other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      done = std::move(other.done);
      new_prime = std::move(other.new_prime);
      check = std::move(other.check);
      use_chinese_remainder = std::move(other.use_chinese_remainder);
      curr_zi_order = std::move(other.curr_zi_order);
      prime_number = std::move(other.prime_number);
      num_eqn = std::move(other.num_eqn);
      n = std::move(other.n);
      type = std::move(other.type);
      zi = std::move(other.zi);
      combined_prime = std::move(other.combined_prime);
    }

    return *this;
  }

  uint BaseReconst::get_num_eqn() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return num_eqn;
  }

  uint BaseReconst::get_prime() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return prime_number;
  }

  //TODO allow for seed with std::srand(std::time(nullptr));
  FFInt BaseReconst::get_rand() {
    return FFInt(pcg32()) + FFInt(1);
  }

  uint BaseReconst::get_zi() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return zi;
  }

  std::vector<uint> BaseReconst::get_zi_order() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end());
  }

  bool BaseReconst::is_done() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return done;
  }

  bool BaseReconst::is_new_prime() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return new_prime;
  }

  mpz_map BaseReconst::convert_to_mpz(const ff_map& coefs) const {
    mpz_map ci_mpz;

    for (const auto & coef : coefs) {
      ci_mpz.insert(std::make_pair(coef.first, mpz_class(coef.second.n)));
    }

    return ci_mpz;
  }

  ff_map BaseReconst::convert_to_ffint(const rn_map& ri) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % FFInt::p);

      if (tmp < 0) tmp = tmp + FFInt::p;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % FFInt::p;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.emplace(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  uint32_t BaseReconst::rotr32(uint32_t x, uint r) {
    return x >> r | x << (-r & 31);
  }

  void BaseReconst::pc32_init(uint64_t seed) {
    state = seed + increment;
    pcg32();
  }

  uint32_t BaseReconst::pcg32() {
    uint64_t x = state;
    uint count = (uint)(x >> 59);

    {
      std::unique_lock<std::mutex> lock_statics(mutex_state);
      state = x * multiplier + increment;
    }

    x ^= x >> 18;
    return rotr32((uint32_t)(x >> 27), count);
  }
}
