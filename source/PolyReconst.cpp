#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"

namespace firefly {

  PolyReconst::PolyReconst (int n_) : n (n_) {
    yi.reserve (5000);
  }

  std::vector<RationalNumber> PolyReconst::reconst() {
    uint64_t first_prime = primes().back();
    combined_prime = first_prime;
    combined_ci = reconst_ff (first_prime);

    for (uint i = primes().size() - 1; i > 0; i--) {
      bool runtest = true;
      for (const auto ci : combined_ci) {
        mpz_class a = ci;
        try {
          guessi.emplace_back (getRationalCoef (a, combined_prime));
        } catch (const std::exception&) {
          runtest = false;
          break;
        }
      }

      uint64_t prime = primes().at (i);
      if (runtest) {
        if (test_guess (prime)) break;
      }
      
      std::vector<mpz_class> ci_tmp = reconst_ff (prime);

      std::pair<mpz_class, mpz_class> p1 (combined_ci.at (0), combined_prime);
      std::pair<mpz_class, mpz_class> p2 (ci_tmp.at (0), prime);

      std::pair<mpz_class, mpz_class> p3 = chineseRemainder (p1, p2);
      combined_ci.at (0) = p3.first;

      for (uint j = 1; j < (uint) combined_ci.size(); j++) {
        p1 = std::make_pair (combined_ci.at (j), combined_prime);
        p2 = std::make_pair (ci_tmp.at (j), prime);

        std::pair<mpz_class, mpz_class> p3j = chineseRemainder (p1, p2);
        combined_ci.at (j) = p3j.first;
      }

      combined_prime = p3.second;
    }

    return guessi;
  }

  std::vector<mpz_class> PolyReconst::reconst_ff (const uint64_t prime) {
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

      if (ai.at (i).n == 0) {
        if (i > breakCondition) {
          bool nonZero = false;

          for (uint j = ai.size(); j > ai.size() - breakCondition; j--) {
            if (ai.at (j - 1).n != 0) {
              nonZero = true;
              break;
            }
          }

          if (!nonZero) break;
        }
      }

      if (i == maxDegree - 1) {
        maxDegree += 5000;
        yi.reserve (maxDegree);
        ai.reserve (maxDegree);
      }
    }

    for (uint i = 0; i < breakCondition; i++) {
      yi.pop_back();
      ai.pop_back();
    }

    yi.reserve (yi.size() + breakCondition);
    return convert_to_mpz (constr_canonical (ai, prime));
  }

  FFInt PolyReconst::comp_ai (const std::vector<FFInt> &ai, int i, int ip, const FFInt &num) {
    if (ip == 0) {
      return num;
    } else {
      if ( (yi.at (i)).n == (yi.at (ip - 1)).n) throw std::runtime_error ("Division by 0 error!");

      return (comp_ai (ai, i, ip - 1, num) - ai.at (ip - 1)) / (yi.at (i) - yi.at (ip - 1));
    }
  }

  std::vector<mpz_class> PolyReconst::convert_to_mpz (const std::vector<FFInt> &ci) const {
    std::vector<mpz_class> ci_mpz;
    ci_mpz.reserve (ci.size());

    for (const auto coef : ci) {
      mpz_class coef_i (coef.n);
      ci_mpz.push_back (coef_i);
    }

    return ci_mpz;
  }

  std::vector<FFInt> PolyReconst::convert_to_ffint (const uint64_t prime) const {
    std::vector<FFInt> gi_ffi;
    gi_ffi.reserve (guessi.size());

    for (const auto gi : guessi) {
      mpz_class tmp (gi.nominator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n (std::stoull (tmp.get_str()), prime);
      tmp = gi.denominator % prime;
      FFInt d (std::stoull (tmp.get_str()), prime);
      gi_ffi.emplace_back (n / d);
    }

    return gi_ffi;
  }


  std::vector<FFInt> PolyReconst::constr_canonical (const std::vector<FFInt> &ai, const uint64_t prime) const {
    if (ai.size() == 0) {
      INFO_MSG ("Polynomial not yet reconstructed or 0.");
      return std::vector<FFInt> {};
    } else if (ai.size() == 1) {
      return ai;
    } else {
      std::vector<FFInt> coef {ai.at (0) };
      Polynomial poly (coef);
      return (poly + iterate_canonical (ai, 1, prime)).coef;
    }
  }

  Polynomial PolyReconst::iterate_canonical (const std::vector<FFInt> &ai, uint i, const uint64_t prime) const {
    if (i < ai.size() - 1) {
      std::vector<FFInt> coef1 { (FFInt (0, prime) - yi.at (i - 1)) *ai.at (i), ai.at (i) };
      std::vector<FFInt> coef2 { FFInt (0, prime) - yi.at (i - 1), FFInt (1, prime) };
      Polynomial poly1 (coef1);
      Polynomial poly2 (coef2);
      return poly1 + poly2 * iterate_canonical (ai, i + 1, prime);
    } else {
      std::vector<FFInt> coef1 { (FFInt (0, prime) - yi.at (i - 1)) *ai.at (i), ai.at (i) };
      Polynomial poly1 (coef1);
      return poly1;
    }
  }

  bool PolyReconst::test_guess (const uint64_t prime) {
    std::vector<FFInt> gi_ffi = convert_to_ffint (prime);
    Polynomial gy (gi_ffi);

    for (uint i = 0; i < std::min(breakCondition, (uint) yi.size()); i++) {
      if (gy.calc (yi.at (i)) != num (prime, yi.at (i))) return false;
    }

    return true;
  }


  FFInt PolyReconst::num (uint64_t p, const FFInt &y) const {
    FFInt a0_0 (3, p);
    FFInt a0_1 (5, p);
    FFInt a1_0 (6, p);
    FFInt a1_1 (7, p);
    FFInt a2 (18, p);
    FFInt a3 (25, p);
    FFInt a4 (30, p);
    FFInt a5 (2, p);
    FFInt a6 (7, p);
    FFInt a7 (100, p);
    FFInt a8 (13, p);
    FFInt exp2 (2, p);
    FFInt exp3 (3, p);
    FFInt exp4 (4, p);
    FFInt exp5 (5, p);
    FFInt exp6 (6, p);
    FFInt exp7 (7, p);
    FFInt exp8 (4000, p);

    return a0_0 / a0_1 - a1_0 / a1_1 * y;
  }
}
