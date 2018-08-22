#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"

namespace firefly {

  PolyReconst::PolyReconst(uint n_) : n(n_) {
    for (uint i = 1; i <= n; i++) {
      std::vector<FFInt> yi;
      yi.reserve(5000);
      yis.insert(std::make_pair(i, std::move(yi)));
      max_deg.insert(std::make_pair(i, -1));
    }
  }

  Polynomial PolyReconst::reconst() {
    uint64_t first_prime = primes().back();
    combined_prime = first_prime;
    const uint zi = n;
    std::vector<uint> chosen_yi(n);
    combined_ci = convert_to_mpz(reconst_ff(zi, first_prime, chosen_yi));

    for (int i = (int) primes().size() - 2; i >= 0; i--) {
      bool runtest = true;

      for (const auto ci : combined_ci) {
        mpz_class a = ci.second;
        try {
          gi.insert(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
        } catch (const std::exception &) {
          runtest = false;
          break;
        }
      }

      uint64_t prime = primes()[i];

      if (runtest) {
        if (test_guess(prime)) break;
      }

      gi.clear();
      yis.clear();

      for (uint i = 1; i <= n; i++) {
        std::vector<FFInt> yi;
        yi.reserve(max_deg[i]);
        yis.insert(std::make_pair(i, std::move(yi)));
      }

      if (i == 0) throw std::runtime_error("Prime numbers not sufficient to reconstruct your coefficients!");

      // use another prime to utilize the Chinese Remainder Theorem to reconstruct the rational
      // coefficients
      mpz_map ci_tmp = convert_to_mpz(reconst_ff(zi, prime, chosen_yi));

      std::pair<mpz_class, mpz_class> p1(combined_ci.begin()->second, combined_prime);
      std::pair<mpz_class, mpz_class> p2(ci_tmp[combined_ci.begin()->first], prime);

      std::pair<mpz_class, mpz_class> p3 = run_chinese_remainder(p1, p2);
      combined_ci.begin()->second = p3.first;

      auto it = combined_ci.begin();

      for (it = ++it; it != combined_ci.end(); ++it) {
        p1 = std::make_pair(it->second, combined_prime);
        p2 = std::make_pair(ci_tmp[it->first], prime);
        p3 = run_chinese_remainder(p1, p2);
        combined_ci[it->first] = p3.first;
      }

      combined_prime = p3.second;
    }

    return Polynomial(gi);
  }

  PolynomialFF PolyReconst::reconst_ff(const uint zi, const uint64_t prime, std::vector<uint> &chosen_yi) {
    std::vector<FFInt> &yi = yis[zi];
    uint maxDegree = max_deg[zi] > 0 ? max_deg[zi] : yi.capacity();
    std::vector<PolynomialFF> ai;
    ai.reserve(maxDegree);

    bool known_prime = yi.capacity() == yi.size();
    if (!known_prime) {
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      yis.insert(std::make_pair(zi, yi));
    }

    if (zi == 1) {
      std::vector<FFInt> chosen_yi_ff;
      for (uint i = 0; i < (uint) chosen_yi.size(); i++) {
        chosen_yi_ff.emplace_back(yis[i + 1][chosen_yi[i]]);
      }

      ai.emplace_back(num(prime, chosen_yi_ff));

    } else {
      ai.emplace_back(reconst_ff(zi - 1, prime, chosen_yi));
    }

    for (uint i = 1; i < maxDegree; i++) {
      if (!known_prime){
        yi.emplace_back(FFInt(std::rand() % prime, prime));
        yis.insert(std::make_pair(zi, yi));
      }

      chosen_yi[zi - 1] = i;

      PolynomialFF fyi;

      bool spuriousPole = true;
      while (spuriousPole) {
        try {
          if (zi == 1) {
            std::vector<FFInt> chosen_yi_ff;
            for (uint j = 0; j < chosen_yi.size(); j++) {
              chosen_yi_ff.emplace_back(yis[j + 1][chosen_yi[j]]);
            }

            fyi = num(prime, chosen_yi_ff);
          } else {
            fyi = reconst_ff(zi - 1, prime, chosen_yi);
          }

          spuriousPole = false;
        } catch (const std::exception &) {
          yi[i] = FFInt(std::rand() % prime, prime);
        }
      }

      spuriousPole = true;

      while (spuriousPole) {
        try {
          ai.emplace_back(comp_ai(zi, ai, fyi, i, i));
          spuriousPole = false;
        } catch (const std::exception &e) {
          yi[i] = FFInt(std::rand() % prime, prime);
        }
      }

      if (ai[i].zero()) {
        if (i > breakCondition) {
          bool nonZero = false;

          for (uint j = ai.size(); j > ai.size() - breakCondition; j--) {
            if (!ai[j - 1].zero()) {
              nonZero = true;
              break;
            }
          }

          if (!nonZero) break;
        }
      }

      if (!known_prime && i == maxDegree) {
        maxDegree += 5000;
        ai.reserve(maxDegree);
      }
    }
    for (uint i = 0; i < breakCondition; i++) {
      ai.pop_back();
    }

    chosen_yi[zi - 1] = 0;

    if(!known_prime) yi.shrink_to_fit();

    if(max_deg[zi] < 0) max_deg[zi] = yi.size();

    return construct_canonical(zi, ai, prime);
  }

  PolynomialFF PolyReconst::comp_ai(const uint zi, const std::vector<PolynomialFF> &ai,
                                    const PolynomialFF &num, int i, int ip) {
    std::vector<FFInt> &yi = yis[zi];

    if (ip == 0) {
      return num;
    } else {
      if (yi[i].n == yi[ip - 1].n) throw std::runtime_error("Division by 0 error!");

      return (comp_ai(zi, ai, num, i, ip - 1) - ai[ip - 1]) / (yi[i] - yi[ip - 1]);
    }
  }

  PolynomialFF PolyReconst::construct_canonical(const uint zi, std::vector<PolynomialFF> &ai,
                                                const uint64_t prime) {
    if (ai.size() == 0) {
      INFO_MSG("Polynomial not yet reconstructed or 0.");
      return PolynomialFF();
    } else if (ai.size() == 1) {
      return ai[0];
    } else {
      return (ai[0] + iterate_canonical(zi, ai, prime, 1));
    }
  }

  PolynomialFF PolyReconst::iterate_canonical(const uint zi,
                                              std::vector<PolynomialFF> &ai,
                                              const uint64_t prime, uint i) {
    std::vector<FFInt> &yi = yis[zi];

    if (i < ai.size() - 1) {
      PolynomialFF poly = ai[i] + iterate_canonical(zi, ai, prime, i + 1);
      return poly.mul(zi) + poly * (FFInt(0, prime) - yi[i - 1]);
    } else {
      return ai[i] * (FFInt(0, prime) - yi[i - 1]) + ai[i].mul(zi);
    }
  }

  bool PolyReconst::test_guess(const uint64_t prime) {
    ff_map gi_ffi = convert_to_ffint(gi, prime);
    PolynomialFF gy(n, gi_ffi);

    std::vector<uint> zero_element(n);

    for (uint i = 0; i < breakCondition; i++) {
      std::vector<FFInt> chosen_yi;

      for (uint j = 0; j < n; j++) {
        chosen_yi.emplace_back(FFInt(std::rand() % prime, prime));
      }

      if (gy.calc(chosen_yi) != num(prime, chosen_yi).coef[zero_element]) return false;
    }

    return true;
  }

  mpz_map PolyReconst::convert_to_mpz(const PolynomialFF &poly) const {
    mpz_map ci_mpz;

    for (const auto &coef : poly.coef) {
      ci_mpz.insert(std::make_pair(coef.first, mpz_class(coef.second.n)));
    }

    return ci_mpz;
  }

  ff_map PolyReconst::convert_to_ffint(const rn_map &ri, const uint64_t prime) const {
    ff_map gi_ffi;

    for (const auto& g_i : ri) {
      mpz_class tmp(g_i.second.numerator % prime);
      mpz_class tmp2(g_i.second.denominator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n(std::stoull(tmp.get_str()), prime);

      tmp = g_i.second.denominator % prime;

      FFInt d(std::stoull(tmp.get_str()), prime);

      gi_ffi.insert(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  PolynomialFF PolyReconst::num(uint64_t prime, const std::vector<FFInt> &chosen_yi) const {
    FFInt y = chosen_yi[0];
    FFInt y2 = chosen_yi[1];
    FFInt y3 = chosen_yi[2];
    FFInt y4 = chosen_yi[3];
    //FFInt y5 = chosen_yi[4];
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
    test = "1234567891098987998798709805302432098098743432098";
    test = test % prime;
    mpz_class ab = test - prime;
    FFInt a7(std::stoull(test.get_str()), prime);
    FFInt a8(13, prime);
    FFInt exp2(2, prime);
    FFInt exp3(3, prime);
    FFInt exp4(4, prime);
    FFInt exp5(5, prime);
    FFInt exp6(6, prime);
    FFInt exp7(7, prime);
    FFInt exp8(500, prime);

    ff_map res;
    res.insert(std::make_pair(std::vector<uint> (n), a2*y2 
    + a3*y3 + a4*y*y2*y3*y4.pow(exp8) + a3/a4*y.pow(exp2)*y3.pow(exp7) 
    - a4*y4.pow(exp2) - a7/a3 * y.pow(exp2)));

    return PolynomialFF(n, res);
  }
}


