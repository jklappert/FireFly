#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <chrono>
#include <algorithm>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::shifted = false;


  RatReconst::RatReconst(uint n_) : n(n_) {
    ti.reserve(300);
    ai.reserve(300);
    combined_prime = FFInt::p;

    if (n > 1) {
      deg_num.emplace_back(-1);
      deg_den.emplace_back(-1);
      curr_zi_order = std::vector<uint> (n, 1);
      curr_zi_order[n - 1] = 0;

      if (!shifted) {
        shift = std::vector<FFInt> (n);

        for (int j = 0; j < n; j++) {
          //shift[j] = FFInt(std::rand() % 100) + FFInt(1);
        }

        shifted = true;
      }
    }
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num) {
    if (!done) {
      // first check if we are done. If not start the reconstruction again using
      // the chinese remainder theorem in combining the previous results
      if (new_prime && n == 1) {
        ti.clear();
        ai.clear();
        ti.emplace_back(new_ti);

        if (rec_rat_coef()) {
          done = test_guess(num);

          if (done) {
            coef_n.clear();
            coef_d.clear();
            combined_di.clear();
            combined_ni.clear();
            combined_prime = 0;
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

      // basic reconstruction algorithm, check if reconstructed function is equal
      // to numeric input and calculate coefficients a_i, check chinese chinese remainder
      // theorem
      zi = 1;

      if (max_deg_num == -1) {
        ti.emplace_back(new_ti);
        const uint i = ti.size() - 1;

        if (i == 0) {
          ai.emplace_back(num);
        } else {
          if (num == comp_fyi(i - 1, i - 1, ti.back())) check = true;

          ai.emplace_back(comp_ai(i, i, num));
        }
      } else {
        uint size = coef_mat.size();

        if (size == 0) {
          coef_mat.reserve(num_eqn);
        }

        // fill matrix
        std::vector<FFInt> eq;
        eq.reserve(num_eqn + 1);

        for (int r = 0; r <= max_deg_num; r++) {
          eq.emplace_back(new_ti.pow(FFInt(r)));
        }

        for (int rp = 1; rp <= max_deg_den; rp++) {
          eq.emplace_back(-new_ti.pow(FFInt(rp)) * num);
        }

        eq.emplace_back(num);

        coef_mat.emplace_back(std::move(eq));

        if (coef_mat.size() == num_eqn) check = true;
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
          std::pair<PolynomialFF, PolynomialFF> canonical = construct_canonical();

          if (max_deg_num == -1) {
            max_deg_num = canonical.first.max_deg()[0];
            max_deg_den = canonical.second.max_deg()[0];
            num_eqn = max_deg_den + max_deg_num + 1;
          }

          std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(canonical);

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

          prime_number ++;
          new_prime = true;
          return;
        } else {
          std::pair<PolynomialFF, PolynomialFF> canonical;

          if (max_deg_num == -1) {
            canonical = construct_canonical();
            PolynomialFF denominator = canonical.second;

            //TODO catch new shift
            if (denominator.min_deg()[0] > 0) {
              INFO_MSG("No constant term in denominator! Trying again with new paramter shift...");

              for (int j = 0; j < n; j++) {
                shift[j] = FFInt(std::rand() % 100) + FFInt(1);
              }

              done = false;
              ai.clear();
              ti.clear();
              return;
            }

            FFInt equializer = FFInt(1) / denominator.coef[denominator.min_deg()];

            canonical.first = canonical.first * equializer;
            canonical.second = denominator * equializer;

            max_deg_num = canonical.first.max_deg()[0];
            max_deg_den = canonical.second.max_deg()[0];
            curr_deg_num = max_deg_num;
            curr_deg_den = max_deg_den;
            num_eqn = max_deg_den + max_deg_num + 1;

            ai.clear();
            ti.clear();
          } else canonical = solve_gauss();

          zi = curr_zi;

          // save the current results to the map to access them later
          for (const auto & coef : canonical.first.coef) {
            const uint deg = coef.first[0];

            if (first_run) {
              PolyReconst rec(n - 1);
              coef_n.emplace(std::make_pair(deg, std::move(rec)));
              deg_num.emplace_back(deg);
              std::vector<uint> zero_deg(n);
              Monomial zero_mon(zero_deg, RationalNumber(0, 1));

              if (deg < max_deg_num) sub_num.emplace(std::make_pair(deg, Polynomial(zero_mon)));
            }

            if (deg <= curr_deg_num) {
              std::vector<uint> key = {deg, zi};
              saved_num_num[curr_zi_order][key] = coef.second;
            }
          }

          for (const auto & coef : canonical.second.coef) {
            const uint deg = coef.first[0];

            if (first_run) {
              PolyReconst rec(n - 1);
              coef_d.emplace(std::make_pair(deg, std::move(rec)));
              deg_den.emplace_back(deg);
              std::vector<uint> zero_deg(n);
              Monomial zero_mon(zero_deg, RationalNumber(0, 1));

              if (deg < max_deg_den) sub_den.emplace(std::make_pair(deg, Polynomial(zero_mon)));
            }

            if (deg <= curr_deg_den) {
              std::vector<uint> key = {deg, zi};
              saved_num_den[curr_zi_order][key] = coef.second;
            }
          }

          if (first_run) {
            std::sort(deg_num.begin(), deg_num.end());
            std::sort(deg_den.begin(), deg_den.end());
            first_run = false;
          }

          // first reconstruct the numerator
          if (curr_deg_num >= 0) {
            PolyReconst rec_num = coef_n[curr_deg_num];
            zi = rec_num.next_zi + 1;
            prime_number = rec_num.prime_number;
            std::vector<uint> tmp_zi_ord = curr_zi_order;

            while (!rec_num.done) {
              try {
                std::vector<uint> key = {(uint) curr_deg_num, zi};
                FFInt food = saved_num_num.at(tmp_zi_ord).at(key);
                // delete unused saved data
                saved_num_num[tmp_zi_ord].erase(key);

                // feed to PolyReconst
                if (curr_deg_num == max_deg_num) rec_num.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
                else {
                  std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
                  yis.emplace(yis.begin(), FFInt(1));
                  FFInt num_subtraction = sub_num[curr_deg_num].convert_to_PolynomialFF().calc(yis);
                  yis.erase(yis.begin());
                  rec_num.feed(yis, food - num_subtraction);
                }

                if (rec_num.next_zi + 1 != zi || n == 2) {
                  zi = rec_num.next_zi + 1;

                  // here one needs to reset everything up the the highest zi value
                  // with ones. The iterator at end has value n - 1 so we caluclate
                  // n - 1 - (n - 1 - zi + 2) = zi - 2 which es equal to the element
                  // before the highest zi value
                  if (prime_number != rec_num.prime_number) {
                    std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1, 1);
                    tmp_zi_ord[n - 1] = rec_num.prime_number;
                    FFInt::p = primes()[rec_num.prime_number];
                  } else {
                    tmp_zi_ord[zi - 2] ++;
                    std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
                  }
                } else tmp_zi_ord[zi - 2] ++;

                prime_number = rec_num.prime_number;
              } catch (std::out_of_range& e) {
                coef_n[curr_deg_num] = rec_num;
                curr_zi = zi;
                curr_zi_order = tmp_zi_ord;
                //std::cout << "catch " << curr_deg_num << " " << zi << " " << tmp_zi_ord[0] << " " << tmp_zi_ord[1] << " " << " " << rec_num.prime_number << " " << prime_number << " " << primes()[prime_number] << " " << FFInt::p << "\n";
                return;
              }

              if (rec_num.done) {
                coef_n[curr_deg_num] = rec_num;
                sub_num.erase(curr_deg_num);
                Polynomial sub_num_pol = rec_num.get_result().homogenize(curr_deg_num).add_shift(shift);
                sub_num_pol -= rec_num.get_result().homogenize(curr_deg_num);

                for (auto & el : sub_num_pol.coefs) {
                  uint tmp_deg = 0;

                  for (auto & n : el.powers) {
                    tmp_deg += n;
                  }

                  sub_num[tmp_deg] += el;
                }

                deg_num.pop_back();
                curr_deg_num = deg_num.back();

                if (curr_deg_num >= 0) {
                  rec_num = coef_n[curr_deg_num];
                  std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1, 1);
                  tmp_zi_ord[n - 1] = 0;
                  zi = rec_num.next_zi + 1;
                  prime_number = rec_num.prime_number;
                  FFInt::p = primes()[prime_number];
                } else break;
              }
            }
          }

          // reconstruct the denominator TODO make the code more elegant
          if (curr_deg_num < 0 && curr_deg_den >= 0) {
            PolyReconst rec_den = coef_d[curr_deg_den];
            zi = rec_den.next_zi + 1;
            prime_number = rec_den.prime_number;
            FFInt::p = primes()[prime_number];
            std::vector<uint> tmp_zi_ord;

            if (first_den_rec) {
              first_den_rec = false;
              tmp_zi_ord = std::vector<uint> (n, 1);
              tmp_zi_ord[n - 1] = 0;
            } else tmp_zi_ord = curr_zi_order;

            while (!rec_den.done) {
              try {
                std::vector<uint> key = {(uint) curr_deg_den, zi};
                FFInt food = saved_num_den.at(tmp_zi_ord).at(key);
                // delete unused saved data
                saved_num_den[tmp_zi_ord].erase(key);

                // feed to PolyReconst
                if (curr_deg_den == max_deg_den) rec_den.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
                else {
                  std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
                  yis.emplace(yis.begin(), FFInt(1));
                  FFInt num_subtraction = sub_den[curr_deg_den].convert_to_PolynomialFF().calc(yis);
                  yis.erase(yis.begin());
                  rec_den.feed(yis, food - num_subtraction);
                }

                if (rec_den.next_zi + 1 != zi || n == 2) {
                  zi = rec_den.next_zi + 1;

                  if (prime_number != rec_den.prime_number) {
                    std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1, 1);
                    tmp_zi_ord[n - 1] = rec_den.prime_number;
                    FFInt::p = primes()[rec_den.prime_number];
                  } else {
                    tmp_zi_ord[zi - 2] ++;
                    std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
                  }
                } else tmp_zi_ord[zi - 2] ++;

                prime_number = rec_den.prime_number;
              } catch (std::out_of_range& e) {
                coef_d[curr_deg_den] = rec_den;
                curr_zi_order = tmp_zi_ord;
                curr_zi = zi;
                //std::cout << "catch den " << curr_deg_den << " " << zi << " " << tmp_zi_ord[0] << " " << tmp_zi_ord[1] << " " << " " << rec_den.prime_number << " " << prime_number << " " << primes()[prime_number] << " " << FFInt::p << "\n";
                return;
              }

              if (rec_den.done) {
                coef_d[curr_deg_den] = rec_den;
                sub_den.erase(curr_deg_den);
                Polynomial sub_den_pol = rec_den.get_result().homogenize(curr_deg_den).add_shift(shift);
                sub_den_pol -= rec_den.get_result().homogenize(curr_deg_den);

                for (auto & el : sub_den_pol.coefs) {
                  uint tmp_deg = 0;

                  for (auto & n : el.powers) {
                    tmp_deg += n;
                  }

                  sub_den[tmp_deg] += el;
                }

                deg_den.pop_back();
                curr_deg_den = deg_den.back();

                if (curr_deg_den >= 0) {
                  rec_den = coef_d[curr_deg_den];
                  std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1, 1);
                  tmp_zi_ord[n - 1] = 0;
                  zi = rec_den.next_zi + 1;
                  prime_number = rec_den.prime_number;
                  FFInt::p = primes()[prime_number];
                } else break;
              }
            }
          }

          done = curr_deg_den == - 1 && curr_deg_num == -1;
          return;
        }
      }
    }
  }

  RationalFunction RatReconst::get_result() {
    if (result.numerator.coefs.empty()) {
      Polynomial numerator;
      Polynomial denominator;

      if (n == 1) {
        numerator = Polynomial(g_ni);
        denominator = Polynomial(g_di);
        g_ni.clear();
        g_di.clear();
      } else {
        for (auto & el : coef_n) {
          numerator += el.second.get_result().homogenize(el.first);
        }

        for (auto & el : coef_d) {
          denominator += el.second.get_result().homogenize(el.first);
        }
      }

      deg_num.clear();
      deg_den.clear();
      curr_zi_order.clear();
      saved_num_num.clear();
      saved_num_den.clear();
      coef_n.clear();
      coef_d.clear();
      numerator.sort();
      denominator.sort();
      result = RationalFunction(numerator, denominator);

      RationalNumber first_coef = result.denominator.coefs[0].coef;

      if (first_coef.numerator != 1 || first_coef.denominator != 1) normalize();
    }

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

  std::pair< PolynomialFF, PolynomialFF > RatReconst::solve_gauss() {
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

      // Make all rows below this one zeroin the current column
      for (uint k = i + 1; k < num_eqn; k++) {
        FFInt c = -coef_mat[k][i] / coef_mat[i][i];

        for (uint j = i; j < num_eqn + 1; j++) {
          if (i == j) coef_mat[k][j] = FFInt(0);
          else coef_mat[k][j] += c * coef_mat[i][j];
        }
      }
    }

    // Solve equation A * x = b for an upper triangular matrix
    std::vector<FFInt> results(num_eqn);

    for (int i = num_eqn - 1; i >= 0; i--) {
      results[i] = coef_mat[i][num_eqn] / coef_mat[i][i];

      for (int k = i - 1; k >= 0; k--) {
        coef_mat[k][num_eqn] -= coef_mat[k][i] * results[i];
      }
    }

    coef_mat.clear();
    // Bring result in canonical form
    ff_map numerator;
    ff_map denominator;

    for (uint i = 0; i <= max_deg_num; i++) {
      if (results[i] != FFInt(0)) {
        std::vector<uint> power = {i};
        numerator.emplace(std::make_pair(std::move(power), results[i]));
      }
    }

    const std::vector<uint> zero_power = {0};
    denominator.emplace(std::make_pair(std::move(zero_power), FFInt(1)));

    for (uint i = 1; i <= max_deg_den; i++) {
      if (results[i + max_deg_num] != FFInt(0)) {
        std::vector<uint> power = {i};
        denominator.emplace(std::make_pair(std::move(power), results[i + max_deg_num]));
      }
    }

    return std::make_pair(PolynomialFF(1, numerator), PolynomialFF(1, denominator));
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
}
