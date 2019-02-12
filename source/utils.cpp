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

  RationalNumber get_rational_coef(const mpz_class& a, const mpz_class& p) {
    mpz_class t = 0;
    mpz_class newt = 1;
    mpz_class tmpt;
    mpz_class r = p;
    mpz_class newr = a;
    mpz_class tmpr;
    mpz_class q;

    while (2 * newr * newr > p) {
      q = r / newr;
      tmpt = t;
      t = newt;
      newt = tmpt - q * newt;
      tmpr = r;
      r = newr;
      newr = tmpr - q * newr;
    }

    if (2 * newt * newt > p || gcd(newr, newt) != 1) throw std::runtime_error("Rational reconstruction failed!");

    if (newt < 0) return RationalNumber(-newr, abs(newt));

    return RationalNumber(newr, newt);
  }

  /* Implementation of MQRR from
   * Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   * by M. Monagan
   */
  RationalNumber get_rational_coef_mqrr(const mpz_class& u, const mpz_class& p) {
//    throw std::runtime_error("Rational reconstruction failed!");
    // set to T so that less than one per mil will be false positive results
    mpz_class T = 16384 * mpz_sizeinbase(p.get_mpz_t(), 2);

    if (u == 0) {
      if (p > T) {
        return RationalNumber(0, 1);
      } else {
        throw std::runtime_error("Rational reconstruction failed!");
      }
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
      throw std::runtime_error("Rational reconstruction failed!");
    }

    if (d < 0) {
      return RationalNumber(-n, abs(d));
    }

    return RationalNumber(n, d);
  }

  std::vector<FFInt> solve_gauss_system(uint num_eqn,
                                        std::vector<std::vector<FFInt>>& coef_mat) {

    // Transform the matrix in upper triangular form
    for (uint i = 0; i < num_eqn; i++) {
      // search for maximum in this column
      FFInt max_el = coef_mat[i][i];
      uint max_row = i;

      for (uint k = i + 1; k < num_eqn; k++) {
        FFInt tmp = coef_mat[k][i];

        if (tmp.n > max_el.n) {
          max_el = tmp;
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (uint k = i; k < num_eqn + 1; k++) {
        FFInt tmp = coef_mat[max_row][k];
        coef_mat[max_row][k] = coef_mat[i][k];
        coef_mat[i][k] = tmp;
      }

      // Make all rows below this one zero in the current column
      for (uint k = i + 1; k < num_eqn; k++) {
        FFInt c = -coef_mat[k][i] / coef_mat[i][i];

        for (uint j = i; j < num_eqn + 1; j++) {
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
    } else{
      ERROR_MSG("Singular system of equations!");
      std::exit(-1);
    }

    return results;
  }

}
