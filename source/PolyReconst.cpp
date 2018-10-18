#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"

namespace firefly {

  PolyReconst::PolyReconst(uint n_) : n(n_) {
    next_zi = 1;
    curr_zi = 1;
    combined_prime = FFInt::p;
    curr_zi_order = std::vector<uint>(n, 1);
  }

  PolyReconst::PolyReconst() {}

  void PolyReconst::feed(const std::vector<FFInt>& new_yis, const FFInt& num) {
    if (!done) {
      // if no yi's/ai's are currently stored, initialize everything
      if (ais.empty() && !use_chinese_remainder) {
        for (uint i = 1; i <= n; i++) {
          std::vector<FFInt> yi;
          std::vector<PolynomialFF> ai;
          yi.reserve(300);
          ai.reserve(300);
          yi.emplace_back(new_yis[i - 1]);
          yis.emplace(std::make_pair(i, std::move(yi)));
          ais.emplace(std::make_pair(i, std::move(ai)));
          max_deg.insert(std::make_pair(i, -1));
        }

        yis[next_zi].pop_back();
      }

      if (new_prime) {
        bool runtest = true;

        for (uint i = 0; i < n; i++) {
          yis[i + 1].clear();
          ais[i + 1].clear();
          yis[i + 1].emplace_back(new_yis[i]);
        }

        for (const auto ci : combined_ci) {
          mpz_class a = ci.second;

          try {
            gi.insert(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
          } catch (const std::exception&) {
            runtest = false;
            break;
          }
        }

        if (runtest) {
          done = test_guess(num);

          if (done) {
            yis.clear();
            ais.clear();
            combined_prime = 0;
            combined_ci.clear();
            max_deg.clear();
            new_prime = false;
            use_chinese_remainder = false;
            return;
          }
        }

        gi.clear();
        next_zi = 1;
        curr_zi = 1;

        if (!use_chinese_remainder) use_chinese_remainder = true;

        new_prime = false;
        yis[next_zi].pop_back();
      }

      //if (yis[next_zi].back() != new_yis[next_zi - 1]) {
      for (uint i = 1; i <= next_zi; i++) {
        yis[i].emplace_back(new_yis[i - 1]);
      }

      //} else throw std::runtime_exception("Division by 0 error!");

      next_zi = 1;
      uint i = yis[next_zi].size() - 1;

      // calc ai's for the lowest stage. If check return otherwise set check
      if (i == 0) {
        std::vector<uint> zero_element(n);
        ff_map zero_map;
        zero_map.emplace(std::make_pair(std::move(zero_element), num));
        ais[next_zi].emplace_back(PolynomialFF(n, zero_map));
      } else {
        std::vector<uint> zero_element(n);
        ff_map zero_map;
        zero_map.emplace(std::make_pair(std::move(zero_element), num));
        ais[next_zi].emplace_back(comp_ai(next_zi, i, i, PolynomialFF(n, zero_map), ais[next_zi]));
      }

      // if the lowest stage ai is zero, combine them into an ai for a higher stage
      // and check if we are done
      if (ais[next_zi][i].zero()) {
        ais[next_zi].pop_back();
        yis[next_zi].pop_back();

        if (curr_zi == 1 && n > 1) curr_zi ++;

        if (curr_zi == 1 && n == 1) check = true;

        if (n > 1) {
          for (uint j = next_zi + 1; j <= curr_zi; j++) {
            if (ais[j].empty()) {
              ais[j].emplace_back(construct_canonical(j - 1, ais[j - 1]));
            } else {
              uint k = yis[j].size() - 1;
              PolynomialFF num_ff = construct_canonical(j - 1, ais[j - 1]);
              ais[j].emplace_back(comp_ai(j, k, k, num_ff, ais[j]));
            }

            ais[j - 1].clear();
            yis[j - 1].clear();

            if (!ais[j].back().zero()) {
              next_zi = j;
              return;
            }

            ais[j].pop_back();
            yis[j].pop_back();

            if (max_deg[j] < 0) max_deg[j] = yis[j].size();

            if (j == curr_zi) {
              if (curr_zi == n) {
                next_zi = j;
                check = true;
                break;
              } else {
                curr_zi ++;
                next_zi = curr_zi;
              }
            }
          }
        }

        if (check && next_zi == curr_zi && curr_zi == n) {
          mpz_map ci_tmp = convert_to_mpz(construct_canonical(n, ais[n]));

          if (!use_chinese_remainder) {
            combined_ci = ci_tmp;
          } else {
            // use another prime to utilize the Chinese Remainder Theorem to reconstruct the rational
            // coefficients

            std::pair<mpz_class, mpz_class> p1;
            std::pair<mpz_class, mpz_class> p2;
            std::pair<mpz_class, mpz_class> p3;

            for (auto it = combined_ci.begin(); it != combined_ci.end(); ++it) {
              p1 = std::make_pair(it->second, combined_prime);
              p2 = std::make_pair(ci_tmp[it->first], FFInt::p);
              p3 = run_chinese_remainder(p1, p2);
              combined_ci[it->first] = p3.first;
            }

            combined_prime = p3.second;
          }

          new_prime = true;
          prime_number ++;
          return;
        }

        return;
      }
    }
  }

  Polynomial PolyReconst::get_result() {
    if (result.coefs.empty()) {
      if (done) {
        result = Polynomial(gi);
        result.sort();
        gi.clear();
      } else {
        PolynomialFF poly = construct_canonical(n, ais[n]);
        rn_map res;

        for (auto & el : poly.coef) {
          res.emplace(std::make_pair(el.first, RationalNumber(el.second.n, 1)));
        }

        return Polynomial(res);
      }
    }

    return result;
  }

  PolynomialFF PolyReconst::comp_ai(const uint zi, int i, int ip,
                                    const PolynomialFF& num, std::vector<PolynomialFF>& ai) {
    std::vector<FFInt>& yi = yis[zi];

    if (ip == 0) return num;

    return (comp_ai(zi, i, ip - 1, num, ai) - ai[ip - 1]) / (yi[i] - yi[ip - 1]);
  }

  PolynomialFF PolyReconst::construct_canonical(const uint zi, std::vector<PolynomialFF>& ai) {
    if (ai.size() == 1) return ai[0];

    return (ai[0] + iterate_canonical(zi, 1, ai));
  }

  PolynomialFF PolyReconst::iterate_canonical(const uint zi, uint i, std::vector<PolynomialFF>& ai) {
    std::vector<FFInt>& yi = yis[zi];

    if (i < ai.size() - 1) {
      PolynomialFF poly = ai[i] + iterate_canonical(zi, i + 1, ai);
      return poly.mul(zi) + poly * (-yi[i - 1]);
    }

    return ai[i] * (-yi[i - 1]) + ai[i].mul(zi);
  }

  bool PolyReconst::test_guess(const FFInt& num) {
    ff_map gi_ffi = convert_to_ffint(gi);
    PolynomialFF gy(n, gi_ffi);
    std::vector<FFInt> chosen_yi(n);

    for (uint i = 1; i <= n; i++) {
      chosen_yi[i - 1] = yis[i][0];
    }

    return gy.calc(chosen_yi) == num;
  }

  mpz_map PolyReconst::convert_to_mpz(const PolynomialFF& poly) const {
    mpz_map ci_mpz;

    for (const auto & coef : poly.coef) {
      ci_mpz.insert(std::make_pair(coef.first, mpz_class(coef.second.n)));
    }

    return ci_mpz;
  }

  ff_map PolyReconst::convert_to_ffint(const rn_map& ri) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % FFInt::p);

      if (tmp < 0) tmp = tmp + FFInt::p;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % FFInt::p;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.insert(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

}
