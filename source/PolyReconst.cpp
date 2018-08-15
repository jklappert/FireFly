#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"

namespace firefly {

  PolyReconst::PolyReconst(int n_) : n(n_) {
    yi.reserve(5000);
  }

  Polynomial PolyReconst::reconst() {
    uint64_t first_prime = primes().back();
    combined_prime = first_prime;
    combined_ci = reconst_ff(first_prime);

    for (int i = (int) primes().size() - 2; i >= 0; i--) {
      bool runtest = true;

      for (const auto ci : combined_ci) {
        mpz_class a = ci;

        try {
          gi.emplace_back(get_rational_coef(a, combined_prime));
        } catch (const std::exception &) {
          runtest = false;
          break;
        }
      }

      uint64_t prime = primes().at(i);

      if (runtest) {
        if (test_guess(prime)) break;
      }

      gi.clear();

      if (i == 0) throw std::runtime_error("Prime numbers not sufficient to reconstruct your coefficients!");

      std::vector<mpz_class> ci_tmp = reconst_ff(prime);

      std::pair<mpz_class, mpz_class> p1(combined_ci.at(0), combined_prime);
      std::pair<mpz_class, mpz_class> p2(ci_tmp.at(0), prime);

      std::pair<mpz_class, mpz_class> p3 = run_chinese_remainder(p1, p2);
      combined_ci.at(0) = p3.first;

      for (uint j = 1; j < (uint) combined_ci.size(); j++) {
        p1 = std::make_pair(combined_ci.at(j), combined_prime);
        p2 = std::make_pair(ci_tmp.at(j), prime);

        std::pair<mpz_class, mpz_class> p3j = run_chinese_remainder(p1, p2);
        combined_ci.at(j) = p3j.first;
      }

      combined_prime = p3.second;
    }

    return Polynomial(gi);
  }

  std::vector<mpz_class> PolyReconst::reconst_ff(const uint64_t prime) {
    uint maxDegree = yi.capacity();
    std::vector<FFInt> ai {};
    ai.reserve(maxDegree + breakCondition);

    yi.clear();
    yi.reserve(maxDegree);

    yi.emplace_back(FFInt(std::rand() % prime, prime));
    ai.emplace_back(num(prime, yi.back()));

    for (uint i = 1; i < maxDegree; i++) {
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      FFInt fyi;

      bool spuriousPole = true;

      while (spuriousPole) {
        try {
          fyi = num(prime, yi.back());
          spuriousPole = false;
        } catch (const std::exception &) {
          yi.pop_back();
          yi.emplace_back(FFInt(std::rand() % prime, prime));
        }
      }

      spuriousPole = true;

      while (spuriousPole) {
        try {
          ai.emplace_back(comp_ai(ai, fyi, i, i));
          spuriousPole = false;
        } catch (const std::exception &) {
          yi.pop_back();
          yi.emplace_back(FFInt(std::rand() % prime, prime));
        }
      }

      if (ai.at(i).n == 0) {
        if (i > breakCondition) {
          bool nonZero = false;

          for (uint j = ai.size(); j > ai.size() - breakCondition; j--) {
            if (ai.at(j - 1).n != 0) {
              nonZero = true;
              break;
            }
          }

          if (!nonZero) break;
        }
      }

      if (i == maxDegree - 1) {
        maxDegree += 5000;
        yi.reserve(maxDegree);
        ai.reserve(maxDegree);
      }
    }

    for (uint i = 0; i < breakCondition; i++) {
      yi.pop_back();
      ai.pop_back();
    }

    yi.reserve(yi.size() + breakCondition);
    return convert_to_mpz(construct_canonical(ai, prime));
  }

  FFInt PolyReconst::comp_ai(const std::vector<FFInt> &ai, const FFInt &num, int i, int ip) {
    if (ip == 0) {
      return num;
    } else {
      if ((yi.at(i)).n == (yi.at(ip - 1)).n) throw std::runtime_error("Division by 0 error!");

      return (comp_ai(ai, num, i, ip - 1) - ai.at(ip - 1)) / (yi.at(i) - yi.at(ip - 1));
    }
  }

  std::vector<FFInt> PolyReconst::construct_canonical(const std::vector<FFInt> &ai, const uint64_t prime) const {
    if (ai.size() == 0) {
      INFO_MSG("Polynomial not yet reconstructed or 0.");
      return std::vector<FFInt> {};
    } else if (ai.size() == 1) {
      return ai;
    } else {
      std::vector<FFInt> coef {ai.at(0) };
      PolynomialFF poly(coef);
      return (poly + iterate_canonical(ai, prime, 1)).coef;
    }
  }

  PolynomialFF PolyReconst::iterate_canonical(const std::vector<FFInt> &ai, const uint64_t prime, uint i) const {
    if (i < ai.size() - 1) {
      std::vector<FFInt> coef1 { (FFInt(0, prime) - yi.at(i - 1)) *ai.at(i), ai.at(i) };
      std::vector<FFInt> coef2 { FFInt(0, prime) - yi.at(i - 1), FFInt(1, prime) };
      PolynomialFF poly1(coef1);
      PolynomialFF poly2(coef2);
      return poly1 + poly2 * iterate_canonical(ai, prime, i + 1);
    } else {
      std::vector<FFInt> coef1 { (FFInt(0, prime) - yi.at(i - 1)) *ai.at(i), ai.at(i) };
      PolynomialFF poly1(coef1);
      return poly1;
    }
  }

  bool PolyReconst::test_guess(const uint64_t prime) {
    std::vector<FFInt> gi_ffi = convert_to_ffint(gi, prime);
    PolynomialFF gy(gi_ffi);

    for (uint i = 0; i < std::min(breakCondition, (uint) yi.size()); i++) {
      if (gy.calc(yi.at(i)) != num(prime, yi.at(i))) return false;
    }

    return true;
  }

  std::vector<mpz_class> PolyReconst::convert_to_mpz(const std::vector<FFInt> &ci) const {
    std::vector<mpz_class> ci_mpz;
    ci_mpz.reserve(ci.size());

    for (const auto coef : ci) {
      mpz_class coef_i(coef.n);
      ci_mpz.emplace_back(coef_i);
    }

    return ci_mpz;
  }

  std::vector<FFInt> PolyReconst::convert_to_ffint(const std::vector<RationalNumber> &ri,
                                                   const uint64_t prime) const {
    std::vector<FFInt> gi_ffi;
    gi_ffi.reserve(ri.size());

    for (const auto g_i : ri) {
      mpz_class tmp(g_i.numerator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n(std::stoull(tmp.get_str()), prime);
      tmp = g_i.denominator % prime;
      FFInt d(std::stoull(tmp.get_str()), prime);
      gi_ffi.emplace_back(n / d);
    }

    return gi_ffi;
  }

  FFInt PolyReconst::num(uint64_t prime, const FFInt &y) const {
    FFInt a0_0(3, prime);
    FFInt a0_1(5, prime);
    FFInt a1_0(6, prime);
    FFInt a1_1(7, prime);
    FFInt a2(18, prime);
    FFInt a3(25, prime);
    FFInt a4(30, prime);
    FFInt a5(2, prime);
    FFInt a6(7, prime);
    mpz_class test;
    test = "1234567891098987984233232323232323232323232323232323232897070743454587987987098053098798708432432098098743432098";
    test = test % prime;
    FFInt a7(std::stoull(test.get_str()), prime);
    FFInt a8(13, prime);
    FFInt exp2(2, prime);
    FFInt exp3(3, prime);
    FFInt exp4(4, prime);
    FFInt exp5(5, prime);
    FFInt exp6(6, prime);
    FFInt exp7(7, prime);
    FFInt exp8(1000, prime);

    return a0_0 + a0_0 * y.pow(exp2);
  }
}
