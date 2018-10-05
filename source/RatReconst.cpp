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

    if (!shifted) {
      shift = std::vector<FFInt> (n);
       if (n > 1) {
         for (auto& el : shift) el = FFInt(std::rand() % 1000000) + FFInt(1);
       }

      shifted = true;
    }

    if (n > 1) {
      deg_num.emplace_back(-1);
      deg_den.emplace_back(-1);
      curr_zi_order = std::vector<uint> (n, 1);
      curr_zi_order[n - 1] = 0;
    }
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord) {
    if (!done) {
      std::vector<uint> tmp_vec;
      std::vector<uint> tmp_vec_rev;
      std::vector<uint> feed_zi_ord_rev;

      if (n > 1) {
        tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
        tmp_vec_rev = tmp_vec;
        feed_zi_ord_rev = feed_zi_ord;
        std::reverse(feed_zi_ord_rev.begin(), feed_zi_ord_rev.end());
        std::reverse(tmp_vec_rev.begin(), tmp_vec_rev.end());
      }

      if (feed_zi_ord_rev >= tmp_vec_rev) {
        if (feed_zi_ord == tmp_vec) {
          // first check if we are done. If not start the reconstruction again using
          // the chinese remainder theorem in combining the previous results
          if (new_prime) {
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
                deg_num.clear();
                deg_den.clear();
                curr_zi_order.clear();
                saved_num_num.clear();
                saved_num_den.clear();
                saved_num_den.clear();
                saved_num_num.clear();
                non_solved_coef_num.clear();
                non_solved_coef_den.clear();
                use_chinese_remainder = false;
                return;
              }
            }

            if (!use_chinese_remainder) use_chinese_remainder = true;

            new_prime = false;
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

            if (size == 0)
              coef_mat.reserve(num_eqn);

            std::vector<FFInt> solved_coef_sub_num {};
            std::vector<FFInt> solved_coef_sub_den {};

            // fill matrix
            std::vector<FFInt> eq;
            eq.reserve(num_eqn + 1);

            std::vector<FFInt> yis;

            if (n > 1)
              yis = std::vector<firefly::FFInt> (curr_zi_order.begin(), curr_zi_order.end() - 1);

            yis.insert(yis.begin(), new_ti);

            for (uint i = 1; i < n; i++) {
              yis[i] = yis[0] * yis[i] + shift[i];
            }

            yis[0] += shift[0];

            std::vector<std::pair<FFInt, FFInt>> t_food = {std::make_pair(new_ti, num)};

            if (n > 1) {
              try {
                std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
                t_food.insert(t_food.end(), saved_ti.at(tmp_vec).begin(), saved_ti.at(tmp_vec).end());
              } catch (std::out_of_range& e) {
                // do nothing
              }
            }

            for (auto food : t_food) {
              FFInt tmp_ti = food.first;
              FFInt tmp_num = food.second;

              for (int r = 0; r <= max_deg_num; r++) {
                if (std::find(non_solved_coef_num.begin(), non_solved_coef_num.end(), r) != non_solved_coef_num.end()) {
                  eq.emplace_back(tmp_ti.pow(FFInt(r)));
                  solved_coef_sub_num.emplace_back(solved_coefs_num[r].convert_to_PolynomialFF().calc(yis));
                } else {
                  FFInt sub = solved_coefs_num[r].convert_to_PolynomialFF().calc(yis);

                  if (sub.n > 0)
                    solved_coef_sub_num.emplace_back(sub);
                }
              }

              for (int rp = min_deg_den + 1; rp <= max_deg_den; rp++) {
                if (std::find(non_solved_coef_den.begin(), non_solved_coef_den.end(), rp) != non_solved_coef_den.end()) {
                  eq.emplace_back(-tmp_ti.pow(FFInt(rp)) * tmp_num);
                  solved_coef_sub_den.emplace_back(solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis));
                } else {
                  FFInt sub = solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis);

                  if (sub.n > 0)
                    solved_coef_sub_den.emplace_back(solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis));
                }
              }

              eq.emplace_back(tmp_ti.pow(FFInt(min_deg_den)) * tmp_num);

              for (auto & solved_coef_num : solved_coef_sub_num) {
                eq.back() += -solved_coef_num;
              }

              for (auto & solved_coef_den : solved_coef_sub_den) {
                eq.back() += solved_coef_den * tmp_num;
              }

              coef_mat.emplace_back(std::move(eq));

              if (coef_mat.size() == num_eqn) {
                check = true;

                if (n > 1) {
                  std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
                  saved_ti.erase(tmp_vec);
                }

                break;
              }
            }
          }

          if (check) {
            check = false;

            // todo not needed anymore. Only if one wants to check twice
            // if (num == comp_fyi(i - 1, i - 1, ti.back())) {

            std::pair<PolynomialFF, PolynomialFF> canonical;

            if (max_deg_num == -1) {
              if (ai.capacity() != ai.size()) {
                ai.shrink_to_fit();
                ti.shrink_to_fit();
              }

              ti.pop_back();
              ai.pop_back();

              canonical = construct_canonical();
              PolynomialFF denominator = canonical.second;

              //TODO catch new shift
              if (n > 1 && denominator.min_deg()[0] > 0) {
                INFO_MSG("No constant term in denominator! Trying again with new paramter shift...");

                for (uint j = 0; j < n; j++) {
                  shift[j] = FFInt(std::rand() % 1000000) + FFInt(1);
                }

                done = false;
                ai.clear();
                ti.clear();
                return;
              }

              max_deg_num = canonical.first.max_deg()[0];
              max_deg_den = canonical.second.max_deg()[0];
              curr_deg_num = max_deg_num;
              curr_deg_den = max_deg_den;
              min_deg_den = canonical.second.min_deg()[0];
              non_solved_coef_num = std::vector<uint> (max_deg_num + 1);
              non_solved_coef_den = std::vector<uint> (max_deg_den);

              FFInt equializer = FFInt(1) / denominator.coef[denominator.min_deg()];

              canonical.first = canonical.first * equializer;
              canonical.second = denominator * equializer;
              std::vector<uint> zero_deg(n);
              Monomial zero_mon(zero_deg, RationalNumber(0, 1));
              solved_coefs_num = std::vector<Polynomial> (max_deg_num + 1, Polynomial(zero_mon));
              solved_coefs_den = std::vector<Polynomial> (max_deg_den - min_deg_den, Polynomial(zero_mon));

              PolynomialFF numerator = canonical.first;
              uint deleted_coefs = 0;
              uint solved_coef_num = 0;
              uint solved_coef_den = 0;

              // check for coefficients which are zero and remove them to save numerical runs
              for (int i = 0; i <= max_deg_num; i++) {
                try {
                  std::vector<uint> pow = {(uint) i};
                  numerator.coef.at(pow);
                  non_solved_coef_num[i - deleted_coefs] = i;
                } catch (std::out_of_range& e) {
                  std::vector<uint> pow(n, 0);
                  pow[0] = i;
                  solved_coefs_num[i] = Polynomial(Monomial(pow, RationalNumber(0, 1)));
                  non_solved_coef_num.erase(non_solved_coef_num.begin() + i - deleted_coefs);
                  deleted_coefs ++;
                  solved_coef_num ++;
                }
              }

              deleted_coefs = 0;

              for (int i = min_deg_den + 1; i <= max_deg_den; i++) {
                try {
                  std::vector<uint> pow = {(uint) i};
                  denominator.coef.at(pow);
                  non_solved_coef_den[i - deleted_coefs - 1 - min_deg_den] = i;
                } catch (std::out_of_range& e) {
                  std::vector<uint> pow(n, 0);
                  pow[0] = i - 1;
                  solved_coefs_den[i - 1] = Polynomial(Monomial(pow, RationalNumber(0, 1)));
                  non_solved_coef_den.erase(non_solved_coef_den.begin() + i - 1 - deleted_coefs);
                  deleted_coefs ++;
                  solved_coef_den ++;
                }
              }

              solved_coefs = solved_coef_num + solved_coef_den;

              num_eqn = max_deg_den + max_deg_num + 1 - min_deg_den - solved_coefs;
              ai.clear();
              ti.clear();
            } else
              canonical = solve_gauss();

            //std::cout << curr_deg_num << " " << canonical.first;

            if (n == 1) {
              std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(canonical);
              combine_primes(tmp);
              prime_number ++;
              saved_ti.clear();
              new_prime = true;
              return;
            } else {
              zi = curr_zi;

              ff_map num_coef = canonical.first.coef;
              ff_map den_coef = canonical.second.coef;

              // save the current results to the map to access them later
              for (uint i = 0; i < non_solved_coef_num.size(); i ++) {
                const uint deg = non_solved_coef_num[i];

                if (first_run) {
                  PolyReconst rec(n - 1);
                  coef_n.emplace(std::make_pair(deg, std::move(rec)));
                  deg_num.emplace_back(deg);

                  if ((int) deg < max_deg_num) {
                    std::vector<uint> zero_deg(n);
                    Monomial zero_mon(zero_deg, RationalNumber(0, 1));
                    sub_num.emplace(std::make_pair(deg, Polynomial(zero_mon)));
                  }
                }

                if ((int) deg <= curr_deg_num) {
                  // this saves some memory since we only need one numerical value
                  // for the constant coefficient
                  if (deg == 0 && first_run) {
                    std::vector<uint> key = {deg, zi};
                    saved_num_num[curr_zi_order][key] = num_coef[ {deg}];
                  } else {
                    std::vector<uint> key = {deg, zi};
                    saved_num_num[curr_zi_order][key] = num_coef[ {deg}];
                  }
                }
              }

              for (uint i = 0; i < non_solved_coef_den.size(); i ++) {
                const uint deg = non_solved_coef_den[i];

                if (first_run) {
                  PolyReconst rec(n - 1);
                  coef_d.emplace(std::make_pair(deg, std::move(rec)));
                  deg_den.emplace_back(deg);

                  if ((int) deg < max_deg_den) {
                    std::vector<uint> zero_deg(n);
                    Monomial zero_mon(zero_deg, RationalNumber(0, 1));
                    sub_den.emplace(std::make_pair(deg, Polynomial(zero_mon)));
                  }

                }

                if ((int) deg <= curr_deg_den) {
                  // this saves some memory since we only need one numerical value
                  // for the constant coefficient
                  if (deg == 0 && first_run) {
                    std::vector<uint> key = {deg, zi};
                    saved_num_den[curr_zi_order][key] = den_coef[ {deg}];
                  } else {
                    std::vector<uint> key = {deg, zi};
                    saved_num_den[curr_zi_order][key] = den_coef[ {deg}];
                  }
                }
              }

              if (first_run) {
                std::sort(deg_num.begin(), deg_num.end());
                std::sort(deg_den.begin(), deg_den.end());
                first_run = false;
              }

              feed_poly();

              return;
            }
          }
        } else if (n > 1 && feed_zi_ord != tmp_vec) {
          try {
            saved_ti.at(feed_zi_ord).emplace_back(std::make_pair(new_ti, num));
          } catch (std::out_of_range& e) {
            std::vector<std::pair<FFInt, FFInt>> tmp_ti = {std::make_pair(new_ti, num)};
            saved_ti[feed_zi_ord] = tmp_ti;
          }
        }
      }
    }
  }

  void RatReconst::feed_poly() {
    // first reconstruct the numerator
    if (curr_deg_num >= 0) {
      PolyReconst rec_num = coef_n[curr_deg_num];
      zi = rec_num.next_zi + 1;

      std::vector<uint> tmp_zi_ord = curr_zi_order;

      while (!rec_num.new_prime) {
        try {
          std::vector<uint> key = {(uint) curr_deg_num, zi};
          FFInt food = saved_num_num.at(tmp_zi_ord).at(key);
          // delete unused saved data
          saved_num_num[tmp_zi_ord].erase(key);

          // feed to PolyReconst
          // since the constant is just a constant, we do not have to get mutliple
          // numerical values to reconstruct the coefficient
          if (curr_deg_num == 0) {
            while (!rec_num.new_prime) {
              if (curr_deg_num == max_deg_num)
                rec_num.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
              else {
                std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
                yis.emplace(yis.begin(), FFInt(1));
                FFInt num_subtraction = sub_num[curr_deg_num].convert_to_PolynomialFF().calc(yis);

                yis.erase(yis.begin());
                rec_num.feed(yis, food - num_subtraction);
              }

              if (rec_num.next_zi + 1 != zi || n == 2) {
                zi = rec_num.next_zi + 1;
                tmp_zi_ord[zi - 2] ++;
                std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
              } else tmp_zi_ord[zi - 2] ++;
            }
          } else {
            if (curr_deg_num == max_deg_num)
              rec_num.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
            else {
              std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
              yis.emplace(yis.begin(), FFInt(1));
              FFInt num_subtraction = sub_num[curr_deg_num].convert_to_PolynomialFF().calc(yis);

              yis.erase(yis.begin());
              rec_num.feed(yis, food - num_subtraction);
            }

            if (rec_num.next_zi + 1 != zi || n == 2) {
              zi = rec_num.next_zi + 1;
              tmp_zi_ord[zi - 2] ++;
              std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
            } else tmp_zi_ord[zi - 2] ++;
          }
        } catch (std::out_of_range& e) {
          coef_n[curr_deg_num] = rec_num;
          curr_zi = zi;
          curr_zi_order = tmp_zi_ord;

          try {
            std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
            std::pair<FFInt, FFInt> key_val = saved_ti.at(tmp_vec).back();
            saved_ti.at(tmp_vec).pop_back();
            feed(key_val.first, key_val.second, tmp_vec);
          } catch (std::out_of_range& e) {
            // do nothing
          }

          return;
        }

        if (rec_num.new_prime) {
          sub_num.erase(curr_deg_num);
          Polynomial sub_num_pol = rec_num.get_result().homogenize(curr_deg_num).add_shift(shift);
          sub_num_pol -= rec_num.get_result().homogenize(curr_deg_num);
          coef_n[curr_deg_num] = rec_num;

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
            tmp_zi_ord[n - 1] = prime_number;
            zi = rec_num.next_zi + 1;
          } else break;
        }
      }
    }

    // reconstruct the denominator TODO make the code more elegant
    if (curr_deg_num < 0 && curr_deg_den >= 0) {
      PolyReconst rec_den = coef_d[curr_deg_den];
      zi = rec_den.next_zi + 1;
      std::vector<uint> tmp_zi_ord;

      if (first_den_rec) {
        first_den_rec = false;
        tmp_zi_ord = std::vector<uint> (n, 1);
        tmp_zi_ord[n - 1] = prime_number;
      } else tmp_zi_ord = curr_zi_order;

      while (!rec_den.new_prime) {
        try {
          std::vector<uint> key = {(uint) curr_deg_den, zi};
          FFInt food = saved_num_den.at(tmp_zi_ord).at(key);
          // delete unused saved data
          saved_num_den[tmp_zi_ord].erase(key);

          // feed to PolyReconst
          // since the constant is just a constant, we do not have to get mutliple
          // numerical values to reconstruct the coefficient
          if (curr_deg_den == 0) {
            while (!rec_den.new_prime) {
              if (curr_deg_den == max_deg_den)
                rec_den.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
              else {
                std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
                yis.emplace(yis.begin(), FFInt(1));
                FFInt num_subtraction = sub_den[curr_deg_den].convert_to_PolynomialFF().calc(yis);

                yis.erase(yis.begin());
                rec_den.feed(yis, food - num_subtraction);
              }

              if (rec_den.next_zi + 1 != zi || n == 2) {
                zi = rec_den.next_zi + 1;
                tmp_zi_ord[zi - 2] ++;
                std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
              } else tmp_zi_ord[zi - 2] ++;
            }
          } else {
            if (curr_deg_den == max_deg_den)
              rec_den.feed(std::vector<FFInt> (tmp_zi_ord.begin(), tmp_zi_ord.end() - 1), food);
            else {
              std::vector<FFInt> yis(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1);
              yis.emplace(yis.begin(), FFInt(1));
              FFInt num_subtraction = sub_den[curr_deg_den].convert_to_PolynomialFF().calc(yis);

              yis.erase(yis.begin());
              rec_den.feed(yis, food - num_subtraction);
            }

            if (rec_den.next_zi + 1 != zi || n == 2) {
              zi = rec_den.next_zi + 1;
              tmp_zi_ord[zi - 2] ++;
              std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - (n + 1 - zi) - 1, 1);
            } else tmp_zi_ord[zi - 2] ++;
          }
        } catch (std::out_of_range& e) {
          coef_d[curr_deg_den] = rec_den;
          curr_zi_order = tmp_zi_ord;
          curr_zi = zi;

          try {
            std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
            std::pair<FFInt, FFInt> key_val = saved_ti.at(tmp_vec).back();
            saved_ti.at(tmp_vec).pop_back();
            feed(key_val.first, key_val.second, tmp_vec);
          } catch (std::out_of_range& e) {
            // do nothing
          }

          return;
        }

        if (rec_den.new_prime) {
          sub_den.erase(curr_deg_den);
          Polynomial sub_den_pol = rec_den.get_result().homogenize(curr_deg_den).add_shift(shift);
          sub_den_pol -= rec_den.get_result().homogenize(curr_deg_den);
          coef_d[curr_deg_den] = rec_den;

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
            tmp_zi_ord[n - 1] = prime_number;
            zi = rec_den.next_zi + 1;
          } else break;
        }
      }
    }

    if (curr_deg_den == - 1 && curr_deg_num == -1) {
      saved_num_num.clear();
      saved_num_den.clear();

      first_run = true;

      // Remove normalization due to the shift
      PolynomialFF numerator;
      PolynomialFF denominator;

      for (auto & el : coef_n) {
        numerator = numerator + el.second.get_result().homogenize(el.first).convert_to_PolynomialFF();
      }

      for (auto & el : coef_d) {
        denominator = denominator + el.second.get_result().homogenize(el.first).convert_to_PolynomialFF();
      }

      coef_n.clear();
      coef_d.clear();

      FFInt first_coef = denominator.coef[denominator.min_deg()];

      // normalize
      FFInt equializer = FFInt(1) / first_coef;

      numerator = numerator * equializer;
      denominator = denominator * equializer;

      std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(std::make_pair(numerator, denominator));

      combine_primes(tmp);

      prime_number ++;
      saved_ti.clear();
      std::fill(curr_zi_order.begin(), curr_zi_order.end() - 1, 1);
      curr_zi_order[n - 1] = prime_number;
      curr_zi = 2;
      zi = 1;
      new_prime = true;
      first_den_rec = true;
    }
  }

  void RatReconst::combine_primes(std::pair<mpz_map, mpz_map>& tmp) {
    if (!use_chinese_remainder) {
      combined_ni = tmp.first;
      combined_di = tmp.second;

      // if the coefficient is not a rational number thus divided by 1,
      // it will not change in the next run and can be omitted to save
      // numerical runs
      if (shift[0].n == 0) {
        mpz_map combined_ni_back = combined_ni;

        for (auto & c_ni : combined_ni_back) {
          try {
            RationalNumber rn = get_rational_coef(c_ni.second, combined_prime);

            if (rn.numerator == c_ni.second && rn.denominator == 1) {
              uint deg = 0;

              for (auto & el : c_ni.first) deg += el;

              remove_ni(deg, c_ni.first, rn);
            }
          } catch (std::exception& e) {
            // do nothing
          }
        }

        mpz_map combined_di_back = combined_di;

        for (auto & c_di : combined_di_back) {
          try {
            RationalNumber rn = get_rational_coef(c_di.second, combined_prime);

            if (rn.numerator == c_di.second && rn.denominator == 1) {
              uint deg = 0;

              for (auto & el : c_di.first) deg += el;

              if ((int) deg != min_deg_den)
                remove_di(deg, c_di.first, rn);
            }
          } catch (std::exception& e) {
            // do nothing
          }
        }

        combined_ni_back.clear();
        combined_di_back.clear();
      }
    } else {
      mpz_map combined_ni_back = combined_ni;
      mpz_map combined_di_back = combined_di;
      mpz_class combined_prime_back = combined_prime;
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

      // Remove already known coefficients from solve algorithm to save numerical runs
      if (shift[0].n == 0) {
        for (auto & c_ni : combined_ni_back) {
          uint deg = 0;

          for (auto & el : c_ni.first) deg += el;

          try {
            RationalNumber last_rn = get_rational_coef(c_ni.second, combined_prime_back);
            RationalNumber curr_rn = get_rational_coef(combined_ni[c_ni.first], combined_prime);

            if (last_rn == curr_rn)
              remove_ni(deg, c_ni.first, curr_rn);
          } catch (std::exception& e) {
            if (c_ni.second == combined_ni[c_ni.first]) {
              RationalNumber rn = RationalNumber(c_ni.second, 1);
              remove_ni(deg, c_ni.first, rn);
            }
          }
        }

        for (auto & c_di : combined_di_back) {
          uint deg = 0;

          for (auto & el : c_di.first) deg += el;

          if ((int) deg != min_deg_den) {
            try {
              RationalNumber last_rn = get_rational_coef(c_di.second, combined_prime_back);
              RationalNumber curr_rn = get_rational_coef(combined_di[c_di.first], combined_prime);

              if (last_rn == curr_rn)
                remove_di(deg, c_di.first, curr_rn);
            } catch (std::exception& e) {
              if (c_di.second == combined_di[c_di.first]) {
                RationalNumber rn = RationalNumber(c_di.second, 1);
                remove_di(deg, c_di.first, rn);
              }
            }
          }
        }

        combined_ni_back.clear();
        combined_di_back.clear();
        combined_prime_back = 0;
      }
    }

    num_eqn = max_deg_den + max_deg_num + 1 - min_deg_den - solved_coefs;
    curr_deg_num = *std::max_element(non_solved_coef_num.begin(), non_solved_coef_num.end());
    curr_deg_den = *std::max_element(non_solved_coef_den.begin(), non_solved_coef_den.end());

    sub_num.clear();
    sub_den.clear();
  }

  RationalFunction RatReconst::get_result() {
    if (result.numerator.coefs.empty()) {
      Polynomial numerator;
      Polynomial denominator;

      numerator = Polynomial(g_ni);
      denominator = Polynomial(g_di);
      g_ni.clear();
      g_di.clear();

      numerator.sort();
      denominator.sort();
      result = RationalFunction(numerator, denominator);

      RationalNumber first_coef = result.denominator.coefs[0].coef;

      if (first_coef.numerator != 1 || first_coef.denominator != 1) result = normalize(result);
    }

    return result;
  }

  bool RatReconst::rec_rat_coef() {
    bool run_test = true;

    for (const auto ci : combined_ni) {
      mpz_class a = ci.second;

      try {
        g_ni[ci.first] = get_rational_coef(a, combined_prime);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    for (const auto ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di[ci.first] = get_rational_coef(a, combined_prime);
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

    const std::vector<uint> min_power = {(uint) min_deg_den};
    denominator.emplace(std::make_pair(std::move(min_power), FFInt(1)));

    uint non_solved_num_size = non_solved_coef_num.size();

    for (uint i = 0; i < non_solved_num_size; i ++) {
      std::vector<uint> power = {non_solved_coef_num[i]};
      numerator.emplace(std::make_pair(std::move(power), results[i]));
    }

    for (uint i = 0; i < non_solved_coef_den.size(); i ++) {
      uint pow = non_solved_coef_den[i];

      if ((int) pow != min_deg_den) {
        std::vector<uint> power = {pow};
        denominator.emplace(std::make_pair(std::move(power), results[i + non_solved_num_size]));
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

  RationalFunction RatReconst::normalize(RationalFunction& rf) {
    RationalNumber equializer = rf.denominator.coefs[0].coef;
    RationalNumber terminator(equializer.denominator, equializer.numerator);

    rf.numerator = rf.numerator * terminator;
    rf.denominator = rf.denominator * terminator;
    return rf;
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
    PolynomialFF g_ny(n, g_ff_ni);
    PolynomialFF g_dy(n, g_ff_di);
    std::vector<FFInt> yis = std::vector<FFInt> (n, FFInt(1));
    yis[0] = ti[0];

    for (uint i = 1; i < n; i++) {
      yis[i] = yis[0] * yis[i] + shift[i];
    }

    yis[0] += shift[0];

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;
  }

  void RatReconst::remove_ni(uint deg, const std::vector<uint>& deg_vec, RationalNumber& rn) {
    g_ni[deg_vec] =  rn;
    combined_ni.erase(deg_vec);
    solved_coefs_num[deg] += Monomial(deg_vec, rn);
    bool remove = true;

    for (auto & c_ni_test : combined_ni) {
      uint deg_test = 0;

      for (auto & el : c_ni_test.first) deg_test += el;

      if (deg_test == deg) remove = false;
    }

    if (remove) {
      solved_coefs ++;
      non_solved_coef_num.erase(std::remove(non_solved_coef_num.begin(),
                                            non_solved_coef_num.end(), deg),
                                non_solved_coef_num.end());
    }
  }

  void RatReconst::remove_di(uint deg, const std::vector<uint>& deg_vec, RationalNumber& rn) {
    g_di[deg_vec] =  rn;
    combined_di.erase(deg_vec);
    solved_coefs_den[deg - (min_deg_den + 1)] += Monomial(deg_vec, rn);
    bool remove = true;

    for (auto & c_di_test : combined_di) {
      uint deg_test = 0;

      for (auto & el : c_di_test.first) deg_test += el;

      if (deg_test == deg) remove = false;
    }

    if (remove) {
      solved_coefs ++;
      non_solved_coef_den.erase(std::remove(non_solved_coef_den.begin(),
                                            non_solved_coef_den.end(), deg),
                                non_solved_coef_den.end());
    }
  }
}
