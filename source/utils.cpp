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
    return std::make_pair(false, RationalNumber(0, 1));
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

    /*for(int i = 0; i < num_eqn; i++){
      for(int j = 0; j <= num_eqn; j++){
        std::cout << coef_mat[i][j] << " ";
      }
      std::cout << "\n";
    }*/

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
      throw std::runtime_error("sing");
    }

    return results;
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
      return a > b;

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

    result.reserve(size - 2);
    std::vector<uint32_t> set = {0, 1};

    for (uint32_t counter = 1; counter < size - 1; ++counter) {
      std::vector<uint32_t> tuple(r);

      uint32_t current_value = counter;

      for (size_t i = 0; i < r; i++) {
        uint32_t digit = current_value % 2;
        tuple[r - i - 1] = set[digit];
        current_value /= 2;
      }

      result.push_back(tuple);
    }

    std::sort(result.begin(), result.end(),
    [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
      return a_grt_b_s(b, a);
    });

    return result;
  }
}
