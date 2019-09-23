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

#include "utils.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"

namespace firefly {
  std::pair<mpz_class, mpz_class> run_chinese_remainder(
    const std::pair<mpz_class, mpz_class>& p1,
    const std::pair<mpz_class, mpz_class>& p2) {
    mpz_class a, n, m1, m2, tmp_c;
    mpz_t tmp;
    n = p1.second * p2.second;
    mpz_init(tmp);
    mpz_invert(tmp, p2.second.get_mpz_t(), p1.second.get_mpz_t());
    tmp_c = mpz_class(tmp);
    m1 = tmp_c * p2.second;
    m2 = (1 - m1) % n;

    if (m2 < 0) m2 = m2 + n;

    a = (m1 * p1.first + m2 * p2.first) % n;

    mpz_clear(tmp);
    return std::pair<mpz_class, mpz_class> (a, n);
  }

  std::pair<bool, RationalNumber> get_rational_coef(const mpz_class& a, const mpz_class& p) {
    mpz_class t = 0;
    mpz_class newt = 1;
    mpz_class tmpt;
    mpz_class r = p;
    mpz_class newr = a;
    mpz_class tmpr;
    mpz_class q;
    bool not_failed = true;

    while (2 * newr * newr > p) {
      q = r / newr;
      tmpt = t;
      t = newt;
      newt = tmpt - q * newt;
      tmpr = r;
      r = newr;
      newr = tmpr - q * newr;
    }

    if (2 * newt * newt > p || gcd(newr, newt) != 1) {
      not_failed = false;
      newr = 1;
      newt = 1;
    }

    if (newt < 0)
      return std::make_pair(not_failed, RationalNumber(-newr, abs(newt)));

    return std::make_pair(not_failed, RationalNumber(newr, newt));
  }

  /* Implementation of MQRR from
   * Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   * by M. Monagan
   */
  std::pair<bool, RationalNumber> get_rational_coef_mqrr(const mpz_class& u, const mpz_class& p) {
    // set to T so that less than one percent will be false positive results
    mpz_class T = 1024 * mpz_sizeinbase(p.get_mpz_t(), 2);
    bool not_failed = true;

    if (u == 0) {
      if (p > T)
        return std::make_pair(not_failed, RationalNumber(0, 1));
      else
        return std::make_pair(false, RationalNumber(0, 1));
    }

    mpz_class n = 0;
    mpz_class d = 0;
    mpz_class t0 = 0;
    mpz_class r0 = p;
    mpz_class t1 = 1;
    mpz_class r1 = u;
    mpz_class tmpr;
    mpz_class tmpt;
    mpz_class q;

    while (r1 != 0 && r0 > T) {
      q = r0 / r1; // mpz_class automatically rounds off

      if (q > T) {
        n = r1;
        d = t1;
        T = q;
      }

      tmpr = r0;
      r0 = r1;
      r1 = tmpr - q * r1;
      tmpt = t0;
      t0 = t1;
      t1 = tmpt - q * t1;
    }

    if (d == 0 || gcd(n, d) != 1) {
      not_failed = false;
      n = 1;
      d = 1;
    }

    if (d < 0)
      return std::make_pair(not_failed, RationalNumber(-n, abs(d)));

    return std::make_pair(not_failed, RationalNumber(n, d));
  }

  bool a_grt_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    for (int i = a.size() - 1; i != -1; --i) {
      if (a[i] == b[i])
        continue;
      else if (a[i] > b[i])
        return true;
      else
        return false;
    }

    return false;
  }

  bool a_grt_b_s(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    uint32_t deg1 = 0;
    uint32_t deg2 = 0;

    for (const auto & el : a) deg1 += el;

    for (const auto & el : b) deg2 += el;

    if (deg1 < deg2)
      return false;
    else if (deg1 == deg2)
      return a_grt_b(b, a);

    return true;
  }

  std::vector<std::vector<uint32_t>> generate_possible_shifts(uint32_t r) {
    std::vector<std::vector<uint32_t>> result;
    uint32_t size = 1;
    uint32_t exp = r;
    uint32_t base = 2;

    while (exp) {
      if (exp & 1)
        size *= base;

      exp >>= 1;

      base *= base;
    }

    result.reserve(size - 1);
    result.emplace_back(std::vector<uint32_t> (r));
    std::vector<uint32_t> set = {0, 1};

    for (uint32_t counter = 1; counter < size - 1; ++counter) {
      std::vector<uint32_t> tuple(r);

      uint32_t current_value = counter;

      for (size_t i = 0; i < r; i++) {
        uint32_t digit = current_value % 2;
        tuple[r - i - 1] = set[digit];
        current_value /= 2;
      }

      result.emplace_back(tuple);
    }

    std::sort(result.begin(), result.end(),
    [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
      return a_grt_b_s(b, a);
    });

    return result;
  }

  PolynomialFF solve_vandermonde_system(std::vector<std::vector<uint32_t>>& degs,
                                        const std::vector<FFInt>& nums,
  const std::vector<FFInt>& val) {
    uint32_t num_eqn = degs.size();
    std::vector<FFInt> result(num_eqn);
    uint32_t n = val.size();

    // calculate base entries of Vandermonde matrix
    std::vector<FFInt> vis;
    vis.reserve(num_eqn);
    std::sort(degs.begin(), degs.end(), std::greater<std::vector<uint32_t>>());

    for (const auto & el : degs) {
      FFInt vi = 1;

      // z_1 is always = 1 which does not matter while determining the coefficient
      for (uint32_t i = 0; i < n; ++i) {
        // curr_zi_ord starts at 1, thus we need to subtract 1 entry
        vi *= val[i].pow(el[i + 1]);
      }

      vis.emplace_back(vi);
    }

    // Initialize the coefficient vector of the master polynomial
    std::vector<FFInt> cis(num_eqn);

    // The coefficients of the master polynomial are found by recursion
    // where we have
    // P(Z) = (Z - v_0)*(Z - v_1)*...*(Z - v_{n-1})
    //      =  c_0 + c_1*Z + ... + Z^n
    cis[num_eqn - 1] = -vis[0];

    for (uint32_t i = 1; i < num_eqn; ++i) {
      for (uint32_t j = num_eqn - 1 - i; j < num_eqn - 1; ++j) {
        cis[j] -= vis[i] * cis[j + 1];
      }

      cis[num_eqn - 1] -= vis[i];
    }

    // Each subfactor in turn is synthetically divided,
    // matrix-multiplied by the right hand-side,
    // and supplied with a denominator (since all vi should be different,
    // there is no additional check if a coefficient in synthetical division
    // leads to a vanishing denominator)
    for (uint32_t i = 0; i < num_eqn; ++i) {
      FFInt t = 1;
      FFInt b = 1;
      FFInt s = nums[num_eqn - 1];

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += nums[j - 1] * b;
        t = vis[i] * t + b;
      }

      result[i] = s / t / vis[i];
    }

    // Bring result in canonical form
    ff_map poly;

    for (uint32_t i = 0; i < num_eqn; ++i) {
      poly.emplace(std::make_pair(degs[i], result[i]));
    }

    return PolynomialFF(n + 1, poly);
  }

#ifdef DEFAULT
  /*uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m) {
    // m must be at most 63 bit.
    uint64_t d = 0;
    uint64_t mp2 = m >> 1;

    for (int i = 0; i != 64; ++i) {
      d = (d > mp2) ? (d << 1) - m : d << 1;

      if (a & 0x8000000000000000ULL) d += b;

      if (d > m) d -= m;

      a <<= 1;
    }

    return d;
  }*/

  uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m) {
    long double x;
    uint64_t c;
    int64_t r;

    if (a >= m) a %= m;

    if (b >= m) b %= m;

    x = a;
    c = x * b / m;
    r = (int64_t)(a * b - c * m) % (int64_t)m;
    return r < 0 ? r + m : r;
  }

  uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t m) {
    // m must be at most 63 bit
    uint64_t res = 1;

    while (exp) {
      if (exp & 1) res = mod_mul(res, base, m);

      base = mod_mul(base, base, m);
      exp >>= 1;
    }

    return res;
  }

  uint64_t mod_inv(uint64_t a, uint64_t m) {
    int64_t t {0};
    int64_t newt {1};
    int64_t tmpt;
    uint64_t r {m};
    uint64_t newr {a};
    uint64_t tmpr;
    uint64_t q;

    while (newr) {
      q = r / newr;
      tmpt = t;
      t = newt;
      newt = tmpt - q * newt;
      tmpr = r;
      r = newr;
      newr = tmpr - q * newr;
    }

    return t < 0 ? t + m : t;
  }
#endif

// Example for the reconstruction of a rational function
  void reconstruct_rational_function() {
    uint32_t n = 3;
    FFInt::set_new_prime(primes()[0]);
    BaseReconst br;
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    br.set_seed(seed);

    RatReconst rec(n);
    rec.set_safe_interpolation();

    // One can set a tag to start from a previously saved run after an interpolation
    // over one prime field was successful
    //rec.set_tag("a");
    // Read in a previously saved file to resume a run from this point
    //rec.start_from_saved_file("ff_save/a_0.txt");

    std::cout << "--------------------------------------------------------------\n";
    std::cout << "Interpolating rational function\n";
    std::cout << "--------------------------------------------------------------\n";
    // Initialize all values. t_yis are the scaled y-values and yis are the
    // unshifted ones
    std::vector<FFInt> t_yis(n);
    t_yis.reserve(n - 1);
    // t is the scaling variable
    FFInt t;

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    // One can use this option to find a sparser shift
    //rec.scan_for_sparsest_shift();

    // Feed loop
    std::vector<FFInt> shift = rec.get_zi_shift_vec();
    /*bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;

    // Generate all possible combinations of shifting variables
    auto shift_vec = generate_possible_shifts(n);

    // Run this loop until a proper shift is found
    while (!found_shift) {
      while (!rec.is_done()) {
        t = rec.get_rand();

        // Add the shift to the scaling variable
        FFInt z1 = t + shift[0];

        for (uint j = 2; j <= n; ++j) {
          t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
        }

        std::vector<FFInt> yis(n);
        yis[0] = z1;

        for (uint j = 1; j < n; ++j) {
          yis[j] = t_yis[j - 1];
        }

        //FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
        //FFInt num = n_eq_1(z1); // example for n = 1
        //FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
        FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
        //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient
        //FFInt num = (FFInt(primes()[0]) + FFInt(primes()[2]) * FFInt(primes()[1]) * yis[0] + 3 * yis[0].pow(2) + yis[1]) / (yis[1]);

        // Feed the algorithm with the current zi_order
        ++count;
        rec.feed(t, num, rec.get_zi_order(), primes_used);
        rec.interpolate();
      }

      found_shift = rec.is_shift_working();
      counter ++;

      if (first) {
        found_shift = false;
        first = false;
        counter = 0;
      }

      rec.set_zi_shift(shift_vec[counter]);
      shift = rec.get_zi_shift_vec();
    }

    rec.set_zi_shift(shift_vec[counter - 1]);
    shift = rec.get_zi_shift_vec();
    rec.accept_shift();
    std::cout << "Total numerical runs to get sparse shift: " << count << ".\n";*/

    // In this loop the whole reconstruction of a function happens
    while (!rec.is_done()) {
      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec.get_prime()) {
        if (!rec.need_shift(primes_used + 1)) {
          rec.disable_shift();
          shift = rec.get_zi_shift_vec();
        }

        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec.get_prime());

        FFInt::set_new_prime(primes()[rec.get_prime()]);
        rec.generate_anchor_points();
        kk = 0;
      }

      // Always set the scaling variable to a random value
      t = rec.get_rand();

      // Add the shift to the scaling variable
      FFInt z1 = t + shift[0];

      for (uint j = 2; j <= n; ++j) {
        t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
      }

      std::vector<FFInt> yis(n);
      yis[0] = z1;

      for (uint j = 1; j < n; ++j) {
        yis[j] = t_yis[j - 1];
      }

      //FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
      //FFInt num = n_eq_1(z1); // example for n = 1
      //FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
      //FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
      //FFInt num = bench_3(yis);
      //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient
      FFInt num = FFInt(primes()[1]) * FFInt(primes()[3]) * (FFInt(primes()[0]) + FFInt(primes()[2]) * FFInt(primes()[1]) * yis[0] + 3 * yis[0].pow(2) + yis[1]) / (FFInt(primes()[1]) + yis[1]);
      //FFInt num = FFInt(primes()[0])*(yis[0]+yis[1]*3)/(yis[2].pow(2)+yis[0]);

      ++kk;
      ++count;

      // Feed the algorithm with the current zi_order
      rec.feed(t, num, rec.get_zi_order(), primes_used);
      rec.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec.get_result();
    std::cout << "--------------------------------------------------------------\n";
  }

  // Example for the reconstruction of a polynomial
  void reconstruct_polynomial() {
    FFInt::set_new_prime(primes()[0]);
    uint32_t n = 2;
    PolyReconst rec_poly(n);
    PolyReconst::set_newton(true);
    PolyReconst::set_bt(true);

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    std::cout << "Interpolating polynomial\n";
    std::cout << "--------------------------------------------------------------\n";
    std::vector<FFInt> yis(n);
    rec_poly.generate_anchor_points();

    // Initialize some counters
    count = 0;
    kk = 0;
    primes_used = 0;

    // Reconstruction loop
    while (!rec_poly.is_done()) {
      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec_poly.get_prime()) {
        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec_poly.get_prime());

        FFInt::set_new_prime(primes()[rec_poly.get_prime()]);
        rec_poly.generate_anchor_points();
        kk = 0;
      }

      yis = rec_poly.get_rand_zi_vec(rec_poly.get_zi_order());

      FFInt num = yis[0] + 3 * yis[0].pow(2) + yis[1];

      ++kk;
      ++count;

      rec_poly.feed(num, rec_poly.get_zi_order(), primes_used);
      rec_poly.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec_poly.get_result();
    std::cout << "--------------------------------------------------------------\n";
  }
}
