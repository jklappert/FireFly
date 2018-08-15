#include "utils.hpp"
#include "Logger.hpp"

namespace firefly {
  /**
   *  zahl und p
   */
  std::pair<mpz_class, mpz_class> run_chinese_remainder(
    const std::pair<mpz_class, mpz_class> &p1,
    const std::pair<mpz_class, mpz_class> &p2) {
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
    //mpz_gcdext(tmp,r,s, p1.second.get_mpz_t(), tmp_c.get_mpz_t());
    //tmp_c = mpz_class(s);
    mpz_clear(tmp);
    return std::pair<mpz_class, mpz_class> (a, n);
  };

  RationalNumber get_rational_coef(const mpz_class &a, const mpz_class &p) {
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

    if (2 * newt * newt > p) throw std::runtime_error("Rational reconstruction failed!");

    if (newt < 0) return RationalNumber(-newr, abs(newt));

    return RationalNumber(newr, newt);
  }

}
