#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"

namespace firefly {

  RatReconst::RatReconst(int n_) : n(n_) {
    //yi.reserve(5000);
  }

  /*RationalFunction RatReconst::reconst() {
    uint64_t first_prime = primes().back();
    combined_prime = first_prime;
    std::pair<std::vector<mpz_class>, std::vector<mpz_class>> tmp  = reconst_ff(first_prime);
    combined_ni = tmp.first;
    combined_di = tmp.second;

    for (int i = (int) primes().size() - 2; i >= 0; i--) {
      bool runtest = true;

      for (const auto ci : combined_ni) {
        mpz_class a = ci;

        try {
          g_ni.emplace_back(get_rational_coef(a, combined_prime));
        } catch (const std::exception &) {
          runtest = false;
          break;
        }
      }

      for (const auto ci : combined_di) {
        mpz_class a = ci;

        try {
          g_di.emplace_back(get_rational_coef(a, combined_prime));
        } catch (const std::exception &) {
          runtest = false;
          break;
        }
      }

      uint64_t prime = primes().at(i);

      if (runtest) {
        if (test_guess(prime)) break;
      }

      if (i == 0) throw std::runtime_error("Prime numbers not sufficient to reconstruct your coefficients!");

      g_ni.clear();
      g_di.clear();
      tmp = reconst_ff(prime);

      // numerator
      std::pair<mpz_class, mpz_class> p1(combined_ni.at(0), combined_prime);
      std::pair<mpz_class, mpz_class> p2(tmp.first.at(0), prime);

      std::pair<mpz_class, mpz_class> p3 = run_chinese_remainder(p1, p2);
      combined_ni.at(0) = p3.first;

      for (uint j = 1; j < (uint) combined_ni.size(); j++) {
        p1 = std::make_pair(combined_ni.at(j), combined_prime);
        p2 = std::make_pair(tmp.first.at(j), prime);

        std::pair<mpz_class, mpz_class> p3j = run_chinese_remainder(p1, p2);
        combined_ni.at(j) = p3j.first;
      }

      // denominator
      for (uint j = 0; j < (uint) combined_di.size(); j++) {
        p1 = std::make_pair(combined_di.at(j), combined_prime);
        p2 = std::make_pair(tmp.second.at(j), prime);

        std::pair<mpz_class, mpz_class> p3j = run_chinese_remainder(p1, p2);
        combined_di.at(j) = p3j.first;
      }

      combined_prime = p3.second;
    }

    return RationalFunction(Polynomial(g_ni), Polynomial(g_di));
  }


  std::pair<std::vector<mpz_class>, std::vector<mpz_class>> RatReconst::reconst_ff(const uint64_t prime)  {
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

      if (fyi == comp_fyi(ai, i - 1, i - 1, yi.back(), prime)) {
        bool nonequal = false;

        for (uint j = 0; j < breakCondition; j++) {
          const FFInt y = FFInt(std::rand() % prime, prime);
          FFInt fy = num(prime, y);

          if (fy != comp_fyi(ai, i - 1, i - 1, y, prime)) {
            nonequal = true;
            break;
          }
        }

        if (!nonequal) {
          yi.pop_back();
          break;
        }
      }

      spuriousPole = true;

      while (spuriousPole) {
        try {
          ai.emplace_back(comp_ai(ai, i, i, fyi));
          spuriousPole = false;
        } catch (const std::exception &) {
          yi.pop_back();
          yi.emplace_back(FFInt(std::rand() % prime, prime));
        }
      }

      if (i == maxDegree - 1) {
        maxDegree += 5000;
        ai.reserve(maxDegree);
        yi.reserve(maxDegree);
      }
    }

    yi.reserve(yi.size() + breakCondition);

    return convert_to_mpz(construct_canonical(ai, prime));
  }

  FFInt RatReconst::comp_ai(std::vector<FFInt> &ai, int i, int ip, const FFInt &num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt aiDum = comp_ai(ai, i, ip - 1, num);

      if (aiDum == ai.at(ip - 1)) throw std::runtime_error("Divide by 0 error!");

      return (yi.at(i) - yi.at(ip - 1)) / (aiDum - ai.at(ip - 1));
    }
  }

  FFInt RatReconst::comp_fyi(std::vector<FFInt> &ai, uint i, uint ip, const FFInt &y, const uint64_t prime) const {
    if (ip == 0) {
      return ai.at(i);
    } else {
      return ai.at(i - ip) + (FFInt(0, prime) - yi.at(i - ip) + y) / comp_fyi(ai, i, ip - 1, y, prime);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::construct_canonical(std::vector<FFInt> &ai, const uint64_t prime) const {
    if (ai.size() == 0) {
      return std::make_pair(PolynomialFF(), PolynomialFF());
    } else if (ai.size() == 1) {
      std::vector<FFInt> coefNom {ai.at(0) };
      std::vector<FFInt> coefDen {FFInt(1, prime) };
      PolynomialFF nom(coefNom);
      PolynomialFF den(coefDen);
      return std::pair<PolynomialFF, PolynomialFF> (nom, den);
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(ai, 1, prime);
      std::vector<FFInt> a0 {ai.at(0) };
      std::vector<FFInt> coefZ {FFInt(0, prime) - yi.at(0), FFInt(1, prime) };
      PolynomialFF constant(a0);
      PolynomialFF zPol(coefZ);
      std::pair<PolynomialFF, PolynomialFF> ratFun(constant * r.first + zPol * r.second,
                                               r.first);
      return normalize(ratFun, prime);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::iterate_canonical(std::vector<FFInt> &ai, uint i, const uint64_t prime) const {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(ai, i + 1 , prime);
      PolynomialFF p1(std::vector<FFInt> {ai.at(i) });
      PolynomialFF p2(std::vector<FFInt> {FFInt(0, prime) - yi.at(i), FFInt(1, prime) });
      return std::pair<PolynomialFF, PolynomialFF> (fnp1.first * p1 + fnp1.second * p2,
                                                fnp1.first);
    } else {
      PolynomialFF p1(std::vector<FFInt> {ai.at(i) });
      PolynomialFF p2(std::vector<FFInt> {FFInt(1, prime) });
      return std::pair<PolynomialFF, PolynomialFF> (p1, p2);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::normalize(std::pair<PolynomialFF, PolynomialFF> &ratFun, const uint64_t prime) const {
    for (auto coef : ratFun.second.coef) {
      if (coef.n != 0) {
        ratFun.first = ratFun.first * (FFInt(1, prime) / coef);
        ratFun.second = ratFun.second * (FFInt(1, prime) / coef);
        return ratFun;
      }
    }

    ERROR_MSG("Could not reconstruct rational function. Still has spurious poles!");
    return ratFun;
  }

  std::pair<std::vector<mpz_class>, std::vector<mpz_class>> RatReconst::convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF> &rf) const {
    std::vector<mpz_class> ci_mpz_1;
    std::vector<mpz_class> ci_mpz_2;
    ci_mpz_1.reserve(rf.first.deg);
    ci_mpz_2.reserve(rf.second.deg);

    for (const auto coef : rf.first.coef) {
      mpz_class coef_i(coef.n);
      ci_mpz_1.emplace_back(coef_i);
    }

    for (const auto coef : rf.second.coef) {
      mpz_class coef_i(coef.n);
      ci_mpz_2.emplace_back(coef_i);
    }

    return std::make_pair(ci_mpz_1, ci_mpz_2);
  }

  std::vector<FFInt> RatReconst::convert_to_ffint(const std::vector<RationalNumber> &ri, const uint64_t prime) const {
    std::vector<FFInt> gi_ffi;
    gi_ffi.reserve(ri.size());

    for (const auto gi : ri) {
      mpz_class tmp(gi.numerator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n(std::stoull(tmp.get_str()), prime);
      tmp = gi.denominator % prime;
      FFInt d(std::stoull(tmp.get_str()), prime);
      gi_ffi.emplace_back(n / d);
    }

    return gi_ffi;
  }

  bool RatReconst::test_guess(const uint64_t prime) {
    std::vector<FFInt> g_ff_ni = convert_to_ffint(g_ni, prime);
    std::vector<FFInt> g_ff_di = convert_to_ffint(g_di, prime);
    PolynomialFF g_ny(g_ff_ni);
    PolynomialFF g_dy(g_ff_di);

    for (uint i = 0; i < std::min(breakCondition, (uint) yi.size()); i++) {
      if ((g_ny.calc(yi.at(i)) / g_dy.calc(yi.at(i))) != num(prime, yi.at(i))) return false;
    }

    return true;
  }

  FFInt RatReconst::num(uint64_t prime, const FFInt &y) const {
    FFInt a0(2, prime);
    FFInt a1(6, prime);
    FFInt a2(1, prime);
    FFInt a3(227, prime);
    FFInt a4(30, prime);
    FFInt a5(2, prime);
    FFInt a6(7, prime);
    mpz_class test;
    test = "1234567891098987984325845845858586879708085484545745874587458787878787874587878787874587878787878798309459864387658765876565987658765765767656565765765765876586586586565865865865865808089897070743454587987987098053098798708432432098098743432098";
    test = test % prime;
    FFInt a7(std::stoull(test.get_str()), prime);
    FFInt a8(13, prime);
    FFInt exp2(2, prime);
    FFInt exp3(3, prime);
    FFInt exp4(4, prime);
    FFInt exp5(5, prime);
    FFInt exp6(6, prime);
    FFInt exp7(7, prime);
    FFInt exp12(12, prime);

    return (a0 / a1 - a3 / a4 * y) / (FFInt(0, prime) - a2 * y.pow(exp12));
  }*/

}
