#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"

namespace firefly {

  RatReconst::RatReconst (int n_) : n (n_) {
    yi.reserve(5000);
  }

  std::pair< std::vector<RationalNumber>, std::vector<RationalNumber>> RatReconst::reconst() {
    uint64_t first_prime = primes().back();
    combined_prime = first_prime;
    std::pair<std::vector<mpz_class>, std::vector<mpz_class>> tmp  = reconst_ff (first_prime);
    combined_ni = tmp.first;
    combined_di = tmp.second;

    for (int i = (int) primes().size() - 2; i >= 0; i--) {
      bool runtest = true;
      for (const auto ci : combined_ni) {
        mpz_class a = ci;
        try {
          g_ni.emplace_back (getRationalCoef (a, combined_prime));
        } catch (const std::exception&) {
          runtest = false;
          break;
        }
      }
      for (const auto ci : combined_di) {
        mpz_class a = ci;
        try {
          g_di.emplace_back (getRationalCoef (a, combined_prime));
        } catch (const std::exception&) {
          runtest = false;
          break;
        }
      }
      uint64_t prime = primes().at (i);
      if (runtest) {
        if (test_guess (prime)) break;
      }

      if(i == 0) throw std::runtime_error("Prime numbers not sufficient to reconstruct your coefficients!");

      g_ni.clear();
      g_di.clear();
      tmp = reconst_ff (prime);

      // numerator
      std::pair<mpz_class, mpz_class> p1 (combined_ni.at (0), combined_prime);
      std::pair<mpz_class, mpz_class> p2 (tmp.first.at (0), prime);

      std::pair<mpz_class, mpz_class> p3 = chineseRemainder (p1, p2);
      combined_ni.at (0) = p3.first;

      for (uint j = 1; j < (uint) combined_ni.size(); j++) {
        p1 = std::make_pair (combined_ni.at (j), combined_prime);
        p2 = std::make_pair (tmp.first.at (j), prime);

        std::pair<mpz_class, mpz_class> p3j = chineseRemainder (p1, p2);
        combined_ni.at (j) = p3j.first;
      }

      // denominator
      for (uint j = 0; j < (uint) combined_di.size(); j++) {
        p1 = std::make_pair (combined_di.at (j), combined_prime);
        p2 = std::make_pair (tmp.second.at (j), prime);

        std::pair<mpz_class, mpz_class> p3j = chineseRemainder (p1, p2);
        combined_di.at (j) = p3j.first;
      }

      combined_prime = p3.second;
    }

    return std::make_pair(g_ni, g_di);
  }


  std::pair<std::vector<mpz_class>, std::vector<mpz_class>> RatReconst::reconst_ff(const uint64_t prime)  {
    uint maxDegree = yi.capacity();
    std::vector<FFInt> ai {};
    ai.reserve (maxDegree + breakCondition);

    yi.clear();
    yi.reserve(maxDegree);
    yi.emplace_back (FFInt (std::rand() % prime, prime));
    ai.emplace_back (num (prime, yi.back()));
    for (uint i = 1; i < maxDegree; i++) {
      yi.emplace_back (FFInt (std::rand() % prime, prime));
      FFInt fyi;
      bool spuriousPole = true;
      while (spuriousPole) {
        try {
          fyi = num (prime, yi.back());
          spuriousPole = false;
        } catch (const std::exception &) {
          yi.pop_back();
          yi.emplace_back (FFInt (std::rand() % prime, prime));
        }
      }
      if (fyi == compFyi (ai, i - 1, yi.back(), prime)) {
        bool nonequal = false;

        for (uint j = 0; j < breakCondition; j++) {
          const FFInt y = FFInt (std::rand() % prime, prime);
          FFInt fy = num (prime, y);

          if (fy != compFyi (ai, i - 1, y, prime)) {
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
          ai.emplace_back (comp_ai (ai, i, i, fyi));
          spuriousPole = false;
        } catch (const std::exception &) {
          yi.pop_back();
          yi.emplace_back (FFInt (std::rand() % prime, prime));
        }
      }

      if (i == maxDegree - 1) {
        maxDegree += 5000;
        ai.reserve(maxDegree);
        yi.reserve(maxDegree);
      }
    }

    yi.reserve (yi.size() + breakCondition);

    return convert_to_mpz (construct_canonical (ai, prime));
  }

  FFInt RatReconst::comp_ai (std::vector<FFInt> &ai, int i, int ip, const FFInt &num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt aiDum = comp_ai (ai, i, ip - 1, num);

      if (aiDum == ai.at (ip - 1)) throw std::runtime_error ("Divide by 0 error!");

      return (yi.at (i) - yi.at (ip - 1)) / (aiDum - ai.at (ip - 1));
    }
  }

  FFInt RatReconst::compFyi (std::vector<FFInt> &ai, int i, const FFInt &y, const uint64_t prime) const {
    return iterateCanonicalNum (ai, i, i, y, prime);
  }

  FFInt RatReconst::iterateCanonicalNum (std::vector<FFInt> &ai, uint i, uint ip, const FFInt &y, const uint64_t prime) const {
    if (ip == 0) {
      return ai.at (i);
    } else {
      return ai.at (i - ip) + (FFInt (0, prime) - yi.at (i - ip) + y) / iterateCanonicalNum (ai, i, ip - 1, y, prime);
    }
  }

  std::pair<Polynomial, Polynomial> RatReconst::construct_canonical(std::vector<FFInt> &ai, const uint64_t prime) const {
    if (ai.size() == 0) {
      return std::make_pair(Polynomial (), Polynomial());
    } else if (ai.size() == 1) {
      std::vector<FFInt> coefNom {ai.at (0) };
      std::vector<FFInt> coefDen {FFInt (1, prime) };
      Polynomial nom (coefNom);
      Polynomial den (coefDen);
      return std::pair<Polynomial, Polynomial> (nom, den);
    } else {
      std::pair<Polynomial, Polynomial> r = iterate_canonical (ai, 1, prime);
      std::vector<FFInt> a0 {ai.at (0) };
      std::vector<FFInt> coefZ {FFInt (0, prime) - yi.at (0), FFInt (1, prime) };
      Polynomial constant (a0);
      Polynomial zPol (coefZ);
      std::pair<Polynomial, Polynomial> ratFun (constant * r.first + zPol * r.second,
                                                r.first);
      return normalize (ratFun, prime);
    }
  }

  std::pair<Polynomial, Polynomial> RatReconst::iterate_canonical (std::vector<FFInt> &ai, uint i, const uint64_t prime) const {
    if (i < ai.size() - 1) {
      std::pair<Polynomial, Polynomial> fnp1 = iterate_canonical (ai, i + 1 , prime);
      Polynomial p1 (std::vector<FFInt> {ai.at (i) });
      Polynomial p2 (std::vector<FFInt> {FFInt (0, prime) - yi.at (i), FFInt (1, prime) });
      return std::pair<Polynomial, Polynomial> (fnp1.first * p1 + fnp1.second * p2,
                                                fnp1.first);
    } else {
      Polynomial p1 (std::vector<FFInt> {ai.at (i) });
      Polynomial p2 (std::vector<FFInt> {FFInt (1, prime) });
      return std::pair<Polynomial, Polynomial> (p1, p2);
    }
  }

  std::pair<Polynomial, Polynomial> RatReconst::normalize (std::pair<Polynomial, Polynomial> &ratFun, const uint64_t prime) const {
    for (auto coef : ratFun.second.coef) {
      if (coef.n != 0) {
        ratFun.first = ratFun.first * (FFInt (1, prime) / coef);
        ratFun.second = ratFun.second * (FFInt (1, prime) / coef);
        return ratFun;
      }
    }

    ERROR_MSG ("Could not reconstruct rational function. Still has spurious poles!");
    return ratFun;
  }

  std::pair<std::vector<mpz_class>, std::vector<mpz_class>> RatReconst::convert_to_mpz (const std::pair<Polynomial, Polynomial> &p) const {
    std::vector<mpz_class> ci_mpz_1;
    std::vector<mpz_class> ci_mpz_2;
    ci_mpz_1.reserve (p.first.deg);
    ci_mpz_2.reserve(p.second.deg);

    for (const auto coef : p.first.coef) {
      mpz_class coef_i (coef.n);
      ci_mpz_1.emplace_back (coef_i);
    }

    for (const auto coef : p.second.coef) {
      mpz_class coef_i (coef.n);
      ci_mpz_2.emplace_back (coef_i);
    }

    return std::make_pair(ci_mpz_1, ci_mpz_2);
  }

  std::vector<FFInt> RatReconst::convert_to_ffint (const uint64_t prime,
                                                   const std::vector<RationalNumber> &guess) const {
    std::vector<FFInt> gi_ffi;
    gi_ffi.reserve (guess.size());

    for (const auto gi : guess) {
      mpz_class tmp (gi.numerator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n (std::stoull (tmp.get_str()), prime);
      tmp = gi.denominator % prime;
      FFInt d (std::stoull (tmp.get_str()), prime);
      gi_ffi.emplace_back (n / d);
    }

    return gi_ffi;
  }

  bool RatReconst::test_guess (const uint64_t prime) {
    std::vector<FFInt> g_ff_ni = convert_to_ffint (prime, g_ni);
    std::vector<FFInt> g_ff_di = convert_to_ffint(prime, g_di);
    Polynomial g_ny (g_ff_ni);
    Polynomial g_dy (g_ff_di);

    for (uint i = 0; i < std::min(breakCondition, (uint) yi.size()); i++) {
      if ((g_ny.calc (yi.at (i))/g_dy.calc (yi.at (i))) != num (prime, yi.at (i))) return false;
    }

    return true;
  }

  FFInt RatReconst::num (uint64_t p, const FFInt &y) const {
    FFInt a0 (2, p);
    FFInt a1 (6, p);
    FFInt a2 (1, p);
    FFInt a3 (227, p);
    FFInt a4 (30, p);
    FFInt a5 (2, p);
    FFInt a6 (7, p);
    mpz_class test;
    test = "1234567891098987984325845845858586879708085484545745874587458787878787874587878787874587878787878798309459864387658765876565987658765765767656565765765765876586586586565865865865865808089897070743454587987987098053098798708432432098098743432098";
    test = test % p;
    FFInt a7 (std::stoull(test.get_str()), p);
    FFInt a8 (13, p);
    FFInt exp2 (2, p);
    FFInt exp3 (3, p);
    FFInt exp4 (4, p);
    FFInt exp5 (5, p);
    FFInt exp6 (6, p);
    FFInt exp7 (7, p);
    FFInt exp12 (100, p);

    return (a0/a1 - a3/a4 * y) / (a2 - a7/a8*y.pow(exp12));
  }

}
