#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <chrono>

namespace firefly {

  RatReconst::RatReconst(uint n_) : n(n_) {
    ti.reserve(5000);
    ai.reserve(5000);
    combined_prime = FFInt::p;
    shift = std::vector<FFInt> (n);
  }

  void RatReconst::feed(const FFInt& new_ti, const std::vector<FFInt>& yis, const FFInt& num) {
    if (!done) {
      // first check if we are done. If not start the reconstruction again using
      // the chinese remainder theorem in combining the previous results
      if (new_prime) {
        if (n == 1) {
          ti.clear();
          ai.clear();
          ti.emplace_back(new_ti);

          if (rec_rat_coef()) {
            done = test_guess(num);

            if (done) {
              new_prime = false;
              use_chinese_remainder = false;
              return;
            }
          }

          g_ni.clear();
          g_di.clear();

          if (!use_chinese_remainder) use_chinese_remainder = true;

          ti.pop_back();
        }

        new_prime = false;
      }

      // basic reconstruction algorithm, check if reconstructed function is equal
      // to numeric input and calculate coefficients a_i, check chinese chinese remainder
      // theorem
      zi = 1;
      ti.emplace_back(new_ti);
      const uint i = ti.size() - 1;

      if (i == 0) {
        ai.emplace_back(num);
      } else {
        if (num == comp_fyi(i - 1, i - 1, ti.back())) check = true;

        ai.emplace_back(comp_ai(i, i, num));
      }

      if (check) {
        check = false;

        // todo not needed anymore. Only if one wants to check twice
        // if (num == comp_fyi(i - 1, i - 1, ti.back())) {
        if (ai.capacity() != ai.size()) {
          ai.shrink_to_fit();
          ti.shrink_to_fit();
        }

        ti.pop_back();
        ai.pop_back();

        if (n == 1) {
          std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(construct_canonical());

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
              p2 = std::make_pair(tmp.first[it->first], FFInt::p);
              p3 = run_chinese_remainder(p1, p2);
              combined_ni[it->first] = p3.first;
            }

            // denominator
            for (auto it = combined_di.begin(); it != combined_di.end(); ++it) {
              p1 = std::make_pair(it->second, combined_prime);
              p2 = std::make_pair(tmp.second[it->first], FFInt::p);
              p3 = run_chinese_remainder(p1, p2);
              combined_di[it->first] = p3.first;
            }

            combined_prime = p3.second;
          }

          new_prime = true;
        } else {
          zi = curr_zi;
          curr_zi = 0;
          new_prime = true;
          done = true;

          std::pair<PolynomialFF, PolynomialFF> canonical = construct_canonical();
          PolynomialFF denominator = canonical.second;

          if (denominator.min_deg()[0] > 0) {
            INFO_MSG("No constant term in denominator! Trying again with new paramter shift...");

            for (int j = 1; j < n; j++) {
              shift[j] = FFInt(std::rand() % 99);
            }

            shifted = true;
            poly_new_prime = false;
            done = false;
            ai.clear();
            ti.clear();
            return;
          }

          FFInt equializer = FFInt(1) / canonical.second.coef[denominator.min_deg()];

          canonical.first = canonical.first * equializer;
          canonical.second = canonical.second * equializer;

          // check whether there is already a map for numerator and denominator
          if (coef_n.empty()) {
            for (const auto coef : canonical.first.coef) {
              PolyReconst rec(n - 1);
              rec.feed(yis, coef.second);
              coef_n.insert(std::make_pair(coef.first[0], std::move(rec)));

              if (curr_zi == 0) curr_zi = rec.next_zi + 1;

              curr_zi = std::min(rec.next_zi + 1, curr_zi);

              if (!rec.new_prime) new_prime = false;

              if (!rec.done) done = false;
            }

            for (const auto coef : canonical.second.coef) {
              PolyReconst rec(n - 1);
              rec.feed(yis, coef.second);
              coef_d.insert(std::make_pair(coef.first[0], std::move(rec)));

              if (curr_zi == 0) curr_zi = rec.next_zi + 1;

              curr_zi = std::min(rec.next_zi + 1, curr_zi);

              if (!rec.new_prime) new_prime = false;

              if (!rec.done) done = false;
            }
          } else {
            for (const auto coef : canonical.first.coef) {
              PolyReconst& rec = coef_n[coef.first[0]];

              if (!rec.done) {
                if (curr_zi == 0) curr_zi = rec.next_zi + 1;

                if (rec.new_prime == poly_new_prime && zi == rec.next_zi + 1) {
                  rec.feed(yis, coef.second);
                }

                curr_zi = std::min(rec.next_zi + 1, curr_zi);
                done = false;

                if (!rec.new_prime) new_prime = false;

              }
            }

            for (const auto coef : canonical.second.coef) {
              PolyReconst& rec = coef_d[coef.first[0]];

              if (!rec.done) {
                if (curr_zi == 0) curr_zi = rec.next_zi + 1;

                if (rec.new_prime == poly_new_prime && zi == rec.next_zi + 1) {
                  rec.feed(yis, coef.second);
                }

                curr_zi = std::min(rec.next_zi + 1, curr_zi);
                done = false;

                if (!rec.new_prime) new_prime = false;
              }
            }
          }

          zi = curr_zi;

          poly_new_prime = new_prime;
          ai.clear();
          ti.clear();
        }

        return;
      }
    }
  }

  RationalFunction RatReconst::get_result() {
    Polynomial numerator;
    Polynomial denominator;

    if (n == 1) {
      numerator = Polynomial(g_ni);
      denominator = Polynomial(g_di);

    } else {
      for (auto & el : coef_n) {
        numerator += el.second.get_result().homogenize(el.first);
      }

      for (auto & el : coef_d) {
        denominator += el.second.get_result().homogenize(el.first);
      }
    }

    numerator.sort();
    denominator.sort();
    result = RationalFunction(numerator, denominator);

    if (shifted) remove_shift();

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
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    for (const auto ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di.emplace(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    return run_test;
  }

  FFInt RatReconst::comp_ai(int i, int ip, const FFInt& num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt ai_i = comp_ai(i, ip - 1, num);

      return (ti[i] - ti[ip - 1]) / (ai_i - ai[ip - 1]);
    }
  }

  FFInt RatReconst::comp_fyi(uint i, uint ip, const FFInt& y) {
    if (ip == 0) {
      return ai[i];
    } else {
      return ai[i - ip] + (-ti[i - ip] + y) / comp_fyi(i, ip - 1, y);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::construct_canonical() {
    if (ai.size() == 1) {
      ff_map numerator_ff;
      std::vector<uint> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[0]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(1);
      FFInt mti = -ti[0];
      std::pair<PolynomialFF, PolynomialFF> ratFun(r.first * ai[0] + r.second * mti + r.second.mul(1),
                                                   r.first);
      return ratFun;
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::iterate_canonical(uint i) {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(i + 1);
      FFInt mti = -ti[i];
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

  std::pair<mpz_map, mpz_map> RatReconst::convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF>& rf) const {
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

  ff_map RatReconst::convert_to_ffint(const rn_map& ri) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % FFInt::p);

      if (tmp < 0) tmp = tmp + FFInt::p;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % FFInt::p;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.emplace(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  bool RatReconst::test_guess(const FFInt& num) {
    ff_map g_ff_ni = convert_to_ffint(g_ni);
    ff_map g_ff_di = convert_to_ffint(g_di);
    PolynomialFF g_ny(1, g_ff_ni);
    PolynomialFF g_dy(1, g_ff_di);
    std::vector<FFInt> yis = {ti[0]};

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;

  }

  void RatReconst::remove_shift() {
    auto start = std::chrono::system_clock::now();
    std::vector<RationalNumber> rn_shift(n);
    std::vector<uint> zero_deg(n);
    std::vector<Polynomial> polys(2);
    polys[0] = result.numerator;
    polys[1] = result.denominator;

    Polynomial tmp_poly;

    for (int i = 0; i < n; i++) {
      rn_shift[i] = RationalNumber(-mpz_class(shift[i].n), 1);
    }

    for (int i = 0; i < 2; i++) {
      for (auto & mon : polys[i].coefs) {
        std::vector<uint> powers = mon.powers;
        Polynomial pow_poly;
        std::vector<uint> decr_power = powers;

        for (int j = 0; j < n; j++) {
          uint deg = powers[j];

          if (deg > 0) {
            std::vector<uint> i_power(n);
            i_power[j] = 1;
            decr_power[j] = 0;
            rn_map sub_shift;
            sub_shift.emplace(std::make_pair(i_power, RationalNumber(1, 1)));
            sub_shift.emplace(std::make_pair(zero_deg, rn_shift[j]));
            Polynomial tmp_pow_poly(sub_shift);
            Polynomial mult_tmp_pow_poly = tmp_pow_poly;

            for (int k = 1; k < deg; k++) {
              tmp_pow_poly = tmp_pow_poly * mult_tmp_pow_poly;
            }

            if (pow_poly.coefs.empty()) pow_poly = tmp_pow_poly;
            else pow_poly = pow_poly * tmp_pow_poly;
          }
        }

        if (!pow_poly.coefs.empty()) tmp_poly += pow_poly * Monomial(decr_power, mon.coef);
      }

      if (i == 0) {
        std::vector<Monomial>& n_coefs = result.numerator.coefs;
        n_coefs.erase(n_coefs.begin() + 1, n_coefs.end());
        result.numerator += tmp_poly;
      }

      if (i == 1) {
        std::vector<Monomial>& d_coefs = result.denominator.coefs;
        d_coefs.erase(d_coefs.begin() + 1, d_coefs.end());
        result.denominator += tmp_poly;
      }

      tmp_poly.clear();
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << diff.count() << "s\n";
  }
}
