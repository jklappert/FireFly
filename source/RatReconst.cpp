#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"

namespace firefly {

  RatReconst::RatReconst(int n_, uint64_t prime) : n(n_) {
    ti.reserve(5000);
    ai.reserve(5000);
    combined_prime = prime;
  }

  void RatReconst::feed(uint64_t prime, const FFInt &new_ti, const std::vector<FFInt> &yis, const FFInt &num) {
    if (!done) {
      // first check if we are done. If not start the reconstruction again using
      // the chinese remainder theorem in combining the previous results
      if (new_prime) {
        ti.clear();
        ai.clear();
        ti.emplace_back(new_ti);

        if (rec_rat_coef()) {
          if (n == 1) {
            done = test_guess(prime, num);

            if (done) return;
          }
        }

        g_ni.clear();
        g_di.clear();

        if (!use_chinese_remainder) use_chinese_remainder = true;

        new_prime = false;
        ti.pop_back();
      }

      // basic reconstruction algorithm, check if reconstructed function is equal
      // to numeric input and calculate coefficients a_i, check chinese chinese remainder
      // theorem
      ti.emplace_back(new_ti);
      const uint i = ti.size() - 1;

      if (check) {
        check = false;

        if (num == comp_fyi(i - 1, i - 1, ti.back(), prime)) {
          if (ai.capacity() != ai.size()) {
            ai.shrink_to_fit();
            ti.shrink_to_fit();
          }

          ti.pop_back();
          ai.pop_back();

          std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(construct_canonical(prime));

          if (!use_chinese_remainder) {
            combined_ni = tmp.first;
            combined_di = tmp.second;
          } else {
            std::pair<mpz_class, mpz_class> p1;
            std::pair<mpz_class, mpz_class> p2;
            std::pair<mpz_class, mpz_class> p3;

            //numerator
            for (auto it = combined_ni.begin(); it != combined_ni.end(); ++it) {
              p1 = std::make_pair(it->second, combined_prime);
              p2 = std::make_pair(tmp.first[it->first], prime);
              p3 = run_chinese_remainder(p1, p2);
              combined_ni[it->first] = p3.first;
            }

            // denominator
            for (auto it = combined_di.begin(); it != combined_di.end(); ++it) {
              p1 = std::make_pair(it->second, combined_prime);
              p2 = std::make_pair(tmp.second[it->first], prime);
              p3 = run_chinese_remainder(p1, p2);
              combined_di[it->first] = p3.first;
            }

            combined_prime = p3.second;
          }

          new_prime = true;
          return;
        }
      }

      if (i == 0) {
        ai.emplace_back(num);
      } else {
        if (num == comp_fyi(i - 1, i - 1, ti.back(), prime)) {
          check = true;
        }

        ai.emplace_back(comp_ai(i, i, num));
      }
    }
  }

  RationalFunction RatReconst::get_result() {
    result = RationalFunction(Polynomial(g_ni), Polynomial(g_di));
    RationalNumber first_coef = result.denominator.coefs[0].coef;

    if (first_coef.numerator != 1 || first_coef.denominator != 1) normalize();

    return result;
  }


  bool RatReconst::rec_rat_coef() {
    bool run_test = true;

    for (const auto ci : combined_ni) {
      mpz_class a = ci.second;

      try {
        g_ni.emplace(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
      } catch (const std::exception &) {
        run_test = false;
        break;
      }
    }

    for (const auto ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di.emplace(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
      } catch (const std::exception &) {
        run_test = false;
        break;
      }
    }

    return run_test;
  }

  FFInt RatReconst::comp_ai(int i, int ip, const FFInt &num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt ai_i = comp_ai(i, ip - 1, num);

      return (ti[i] - ti[ip - 1]) / (ai_i - ai[ip - 1]);
    }
  }

  FFInt RatReconst::comp_fyi(uint i, uint ip, const FFInt &y, const uint64_t prime) {
    if (ip == 0) {
      return ai[i];
    } else {
      return ai[i - ip] + (FFInt(0) - ti[i - ip] + y) / comp_fyi(i, ip - 1, y, prime);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::construct_canonical(const uint64_t prime) const {
    if (ai.size() == 1) {
      ff_map numerator_ff;
      std::vector<uint> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[0]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(1, prime);
      FFInt mti = FFInt(0) - ti[0];
      std::pair<PolynomialFF, PolynomialFF> ratFun(r.first * ai[0] + r.second * mti + r.second.mul(1),
                                                   r.first);
      return ratFun;
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::iterate_canonical(uint i, const uint64_t prime) const {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(i + 1 , prime);
      FFInt mti = FFInt(0) - ti[i];
      return std::pair<PolynomialFF, PolynomialFF> (fnp1.first * ai[i] + fnp1.second.mul(1) + fnp1.second * mti,
                                                    fnp1.first);
    } else {
      ff_map numerator_ff;
      std::vector<uint> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[i]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    }
  }

  void RatReconst::normalize() {
    RationalNumber equializer = result.denominator.coefs[0].coef;
    RationalNumber terminator(equializer.denominator, equializer.numerator);

    result.numerator = result.numerator * terminator;
    result.denominator = result.denominator * terminator;
  }

  std::pair<mpz_map, mpz_map> RatReconst::convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF> &rf) const {
    mpz_map ci_mpz_1;
    mpz_map ci_mpz_2;

    for (const auto coef : rf.first.coef) {
      ci_mpz_1.emplace(coef.first, mpz_class(coef.second.n));
    }

    for (const auto coef : rf.second.coef) {
      ci_mpz_2.emplace(coef.first, mpz_class(coef.second.n));
    }

    return std::make_pair(ci_mpz_1, ci_mpz_2);
  }

  ff_map RatReconst::convert_to_ffint(const rn_map &ri, const uint64_t prime) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % prime);

      if (tmp < 0) tmp = tmp + prime;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % prime;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.insert(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  bool RatReconst::test_guess(const uint64_t prime, const FFInt &num) {
    ff_map g_ff_ni = convert_to_ffint(g_ni, prime);
    ff_map g_ff_di = convert_to_ffint(g_di, prime);
    PolynomialFF g_ny(1, g_ff_ni);
    PolynomialFF g_dy(1, g_ff_di);
    std::vector<FFInt> yis = {ti[0]};

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;

  }

}
