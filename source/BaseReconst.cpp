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

#include "BaseReconst.hpp"
#include "utils.hpp"

namespace firefly {
  uint64_t const BaseReconst::multiplier = 6364136223846793005u;
  uint64_t const BaseReconst::increment  = 1442695040888963407u;
  uint64_t BaseReconst::state = 0x4d595df4d0f33173;
  uint64_t BaseReconst::splitmix64_state = 0x4d595df4d0f33173;
  uint64_t BaseReconst::s[4];
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

  uint32_t BaseReconst::get_num_eqn() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return num_eqn;
  }

  uint32_t BaseReconst::get_prime() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return prime_number;
  }

  FFInt BaseReconst::get_rand() {
    FFInt rand(xoshiro256ss()/*pcg32()*/);

    if (rand != 0)
      return rand;
    else {
      while (rand == 0) {
        rand = FFInt(xoshiro256ss()/*pcg32()*/);
      }

      return rand;
    }
  }

  FFInt BaseReconst::get_rand_64() {
    FFInt rand(xoshiro256ss()/*pcg32()*/);

    if (rand != 0)
      return rand;
    else {
      while (rand == 0) {
        rand = FFInt(xoshiro256ss()/*pcg32()*/);
      }

      return rand;
    }
  }

  void BaseReconst::set_seed(uint64_t seed) {
    pc32_init(seed);
    xoshiro256ss_init(seed);
  }

  uint32_t BaseReconst::get_zi() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return zi;
  }

  std::vector<uint32_t> BaseReconst::get_zi_order() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return std::vector<uint32_t>(curr_zi_order.begin(), curr_zi_order.end());
  }

  bool BaseReconst::is_done() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return done;
  }

  std::pair<bool, uint32_t> BaseReconst::get_done_and_prime() const {
    std::unique_lock<std::mutex> lock(mutex_status);
    return std::make_pair(done, prime_number);
  }

  bool BaseReconst::is_new_prime() const {
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
      FFInt n(g_i.second.numerator);
      FFInt d(g_i.second.denominator);
      gi_ffi.emplace(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  uint32_t BaseReconst::rotr32(uint32_t x, uint32_t r) {
    return x >> r | x << (-r & 31);
  }

  void BaseReconst::pc32_init(uint64_t seed) {
    {
      std::unique_lock<std::mutex> lock_statics(mutex_state);
      state = seed + increment;
    }
    pcg32();
  }

  uint32_t BaseReconst::pcg32() {
    uint64_t x;
    uint32_t count;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_state);
      x = state;
      count = static_cast<uint32_t>((x >> 59));
      state = x * multiplier + increment;
    }

    x ^= x >> 18;
    return rotr32(static_cast<uint32_t>((x >> 27)), count);
  }

  uint64_t BaseReconst::rol64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }

  uint64_t BaseReconst::xoshiro256ss() {
    uint64_t result;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_state);
      result = rol64(s[1] * 5, 7) * 9;
      const uint64_t t = s[1] << 17;

      s[2] ^= s[0];
      s[3] ^= s[1];
      s[1] ^= s[2];
      s[0] ^= s[3];

      s[2] ^= t;
      s[3] = rol64(s[3], 45);

    }
    return result;
  }

  // only used for seeding xoshiro256ss
  uint64_t BaseReconst::splitmix64() {
    uint64_t result = splitmix64_state;

    splitmix64_state = result + 0x9E3779B97f4A7C15;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
  }

  void BaseReconst::xoshiro256ss_init(uint64_t seed) {
    {
      std::unique_lock<std::mutex> lock_statics(mutex_state);
      splitmix64_state = seed + increment;
      s[0] = splitmix64();
      s[1] = splitmix64();
      s[2] = splitmix64();
      s[3] = splitmix64();
    }
  }
}
