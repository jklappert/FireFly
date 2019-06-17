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
#include "Logger.hpp"

namespace firefly {
  /**
   *  zahl und p
   */
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

  std::vector<FFInt> solve_gauss_system(uint32_t num_eqn,
  std::vector<std::vector<FFInt>>& coef_mat) {

    // Transform the matrix in upper triangular form
    for (uint32_t i = 0; i < num_eqn; ++i) {
      // search for maximum in this column
      FFInt max_el = coef_mat[i][i];
      uint32_t max_row = i;

      for (uint32_t k = i + 1; k < num_eqn; ++k) {
        FFInt tmp = coef_mat[k][i];

        if (tmp.n > max_el.n) {
          max_el = tmp;
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (uint32_t k = i; k < num_eqn + 1; ++k) {
        FFInt tmp = coef_mat[max_row][k];
        coef_mat[max_row][k] = coef_mat[i][k];
        coef_mat[i][k] = tmp;
      }

      // Make all rows below this one zero in the current column
      for (uint32_t k = i + 1; k < num_eqn; ++k) {
        FFInt c = -coef_mat[k][i] / coef_mat[i][i];

        for (uint32_t j = i; j < num_eqn + 1; ++j) {
          if (i == j) coef_mat[k][j] = FFInt(0);
          else coef_mat[k][j] += c * coef_mat[i][j];
        }
      }
    }

    std::vector<FFInt> results(num_eqn);

    if (coef_mat[num_eqn - 1][num_eqn - 1] != 0) {
      // Solve equation A * x = b for an upper triangular matrix
      for (int i = num_eqn - 1; i >= 0; i--) {
        results[i] = coef_mat[i][num_eqn] / coef_mat[i][i];

        for (int k = i - 1; k >= 0; k--) {
          coef_mat[k][num_eqn] -= coef_mat[k][i] * results[i];
        }
      }
    } else {
      ERROR_MSG("Singular system of equations!");
      throw std::runtime_error("Singular system of equations!");
    }

    return results;
  }

  PolynomialFF solve_vandermonde_system(std::vector<std::vector<uint32_t>>& degs,
                                        const std::vector<std::pair<FFInt, uint32_t>>& nums,
  const std::vector<FFInt> val) {
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
      FFInt s = nums[num_eqn - 1].first;

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += nums[j - 1].first * b;
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

  uint64_t xorshift64star() {
    uint64_t x = xorshift64star_state;
    x ^= x >> 12; // a
    x ^= x << 25; // b
    x ^= x >> 27; // c
    xorshift64star_state = x;
    return x * 0x2545F4914F6CDD1D;
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
}
