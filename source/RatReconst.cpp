//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

#include "Logger.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include "DenseSolver.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sys/stat.h>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::need_prime_shift = false;
  bool RatReconst::set_singular_system = false;
  ff_pair_map RatReconst::rand_zi {};
  std::mutex RatReconst::mutex_statics;
  std::vector<uint32_t> RatReconst::curr_shift {};

  RatReconst::RatReconst(uint32_t n_) {
    n = n_;
    type = RAT;
    const_den = 0;
    std::unique_lock<std::mutex> lock_status(mutex_status);

    combined_prime = FFInt::p;

    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (shift.empty()) {
        shift = std::vector<FFInt> (n);

        if (n > 1) {
          for (auto & el : shift) el = FFInt(xorshift64star());

          curr_shift = std::vector<uint32_t> (n, 1);
        }
      }
    }

    t_interpolator = ThieleInterpolator();

    if (n > 1) {
      curr_zi_order_num = std::vector<uint32_t> (n - 1, 1);
      curr_zi_order_den = std::vector<uint32_t> (n - 1, 1);
      curr_zi_order = std::vector<uint32_t> (n - 1, 1);
      lock_status.unlock();
      // add a zero to both polynomials to do arithmetics
      ff_map zero_deg {};
      zero_deg.emplace(std::make_pair(std::vector<uint32_t> (n), 0));
      solved_num = PolynomialFF(n, zero_deg);
      solved_den = PolynomialFF(n, zero_deg);

      // fill in the rand_vars for zi_order = 1
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (rand_zi.empty()) {
        lock_statics.unlock();
        generate_anchor_points();
      }
    }
  }

  void RatReconst::scan_for_sparsest_shift() {
    scan = true;

    if (n == 1) {
      ERROR_MSG("You should never want to shift a univariate rational function.");
      std::exit(-1);
    }
  }

  void RatReconst::set_zi_shift(const std::vector<uint32_t>& shifted_zis) {
    {
      std::unique_lock<std::mutex> lock(mutex_statics);
      shift = std::vector<FFInt> (n);
      curr_shift = shifted_zis;

      for (uint32_t i = 0; i < n; ++i) {
        if (shifted_zis[i] == 1)
          shift[i] = FFInt(xorshift64star());
      }
    }
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (!done && fed_prime == prime_number)
      queue.emplace(std::make_tuple(new_ti, num, feed_zi_ord));
  }

  bool RatReconst::interpolate() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (is_interpolating || queue.empty()) return true;
    else {
      is_interpolating = true;

      while (!queue.empty()) {
        auto food = queue.front();
        queue.pop();
        lock.unlock();
        interpolate(std::get<0>(food), std::get<1>(food), std::get<2>(food));

        while (saved_ti.find(curr_zi_order) != saved_ti.end()) {
          /*
          * If not finished, check if we can use some saved runs
          */
          if (saved_ti.find(curr_zi_order) != saved_ti.end()) {
            std::pair<FFInt, FFInt> key_val = saved_ti[curr_zi_order].back();
            saved_ti[curr_zi_order].pop_back();
            interpolate(key_val.first, key_val.second, curr_zi_order);
          }
        }

        lock.lock();
      }
    }

    is_interpolating = false;
    return false;
  }

  void RatReconst::interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord) {
    if (!done) {

      // Compare if the food is the expected food; if not, store it for later use
      if (feed_zi_ord == curr_zi_order) {
        // first check if we are done. If not start the interpolation again using
        // the chinese remainder theorem in combining the previous results
        if (new_prime)
          if (check_if_done(num, new_ti))
            return;

        // basic interpolation algorithm, check if interpolation function is equal
        // to numeric input and calculate coefficients a_i, check chinese chinese remainder
        // theorem
        {
          std::unique_lock<std::mutex> lock(mutex_status);

          if (prime_number < interpolations) zi = 1;
        }

        if (max_deg_num == -1)
          check = t_interpolator.add_point(num, new_ti);
        else {
          std::vector<std::pair<FFInt, FFInt>> t_food = {std::make_pair(new_ti, num)};

          // Prepare food for Gauss system
          if (n > 1) {
            if (saved_ti.find(curr_zi_order) != saved_ti.end()) {
              t_food.insert(t_food.end(), saved_ti[curr_zi_order].begin(), saved_ti[curr_zi_order].end());
              saved_ti.erase(curr_zi_order);
            }
          }

          // Iterate through all feeds and build the uni/multivariate Gauss
          // system
          for (auto food : t_food) {
            FFInt tmp_ti = food.first;
            FFInt tmp_num = food.second;

            // Get yi's for the current feed
            std::vector<FFInt> yis;

            if (n > 1)
              yis = get_rand_zi_vec(curr_zi_order);

            yis.emplace(yis.begin(), 1);

            // build Gauss system for univariate reconstruction needed for
            // multivariate rational functions
            if (prime_number < interpolations)
              build_uni_gauss(tmp_ti, tmp_num, yis);
            else
              build_homogenized_multi_gauss(tmp_ti, tmp_num, yis);

            if (num_eqn == coef_mat.size()) {
              check = true;
              break;
            }
          }
        }

        if (check) {
          check = false;

          std::pair<ff_map, ff_map> canonical;

          // If the maximal/minimal degree of the polynomials are not set
          // determine them and save all information
          if (max_deg_num == -1) {

            canonical = t_interpolator.get_result();
            PolynomialFF numerator = PolynomialFF(1, canonical.first);
            PolynomialFF denominator = PolynomialFF(1, canonical.second);

            if (scan) {
              {
                std::unique_lock<std::mutex> lock(mutex_status);
                done = true;
              }

              uint32_t tmp_max_deg_num = numerator.max_deg()[0];
              uint32_t tmp_max_deg_den = denominator.max_deg()[0];

              {
                std::unique_lock<std::mutex> lock_status(mutex_status);
                shift_works = false;
              }

              if (curr_shift == std::vector<uint32_t> (n, 1)) {
                all_shift_max_degs = {tmp_max_deg_num, tmp_max_deg_den};
              } else if (tmp_max_deg_num == all_shift_max_degs[0] && tmp_max_deg_den == all_shift_max_degs[1]) {
                if (denominator.min_deg()[0] == 0) {
                  normalize_to_den = true;
                  start_deg_den = 1;
                  start_deg_num = 0;
                } else {
                  normalize_to_den = false;
                  start_deg_num = 1;
                  start_deg_den = 0;
                }

                shift_works = true;
              }

              t_interpolator = ThieleInterpolator();
              return;
            }

            max_deg_num = numerator.max_deg()[0];
            max_deg_den = denominator.max_deg()[0];

            if (n == 1)
              normalizer_deg = denominator.min_deg();

            if (normalize_to_den) {
              curr_deg_num = max_deg_num;

              if (max_deg_den > 0)
                curr_deg_den = max_deg_den;
            } else {
              curr_deg_den = max_deg_den;

              if (max_deg_num > 0)
                curr_deg_num = max_deg_num;
            }

            FFInt equalizer;

            if (normalize_to_den)
              equalizer = FFInt(1) / denominator.coefs[denominator.min_deg()];
            else
              equalizer = FFInt(1) / numerator.coefs[numerator.min_deg()];

            canonical.first = (numerator * equalizer).coefs;
            canonical.second = (denominator * equalizer).coefs;

            for (uint32_t i = start_deg_num; i < (uint32_t) max_deg_num; ++i) {
              if (canonical.first.find( {i}) == canonical.first.end())
                tmp_solved_coefs_num ++;
              else if (i == 0)
                tmp_sol_const_num = 1;
            }

            for (uint32_t i = start_deg_den; i < (uint32_t) max_deg_den; ++i) {
              if (canonical.second.find( {i}) == canonical.second.end())
                tmp_solved_coefs_den ++;
              else if (i == 0)
                tmp_sol_const_den = 1;
            }

            // set number of equations needed for univariate rational function
            // reconstruction needed for multivariate polynomial feed
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              num_eqn = max_deg_den + max_deg_num + 1
                        - tmp_solved_coefs_num - tmp_solved_coefs_den
                        - tmp_sol_const_num - tmp_sol_const_den;
            }

            t_interpolator = ThieleInterpolator();
          } else if (prime_number < interpolations)
            canonical = solve_gauss();
          else
            canonical = solve_homogenized_multi_gauss();

          if (n == 1) {
            combine_primes(canonical.first, canonical.second);
            return;
          } else if (prime_number < interpolations) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi = curr_zi;
            }

            // remove zero degrees
            if (first_run) {
              ff_map num_coef = canonical.first;
              ff_map den_coef = canonical.second;

              if (num_coef.size() > 1) {
                for (const auto & el : num_coef) {
                  if (el.second == 0) {
                    canonical.first.erase(el.first);
                    tmp_solved_coefs_num ++;
                  }
                }
              }

              for (const auto & el : den_coef) {
                if (el.second == 0) {
                  canonical.second.erase(el.first);
                  tmp_solved_coefs_den ++;
                }
              }
            }

            // save the current results to the map to access them later
            for (const auto & el : canonical.first) {
              uint32_t deg = el.first[0];

              if (first_run) {
                if (!normalize_to_den && deg == 0)
                  continue;

                PolyReconst rec(n - 1, deg, true);
                coef_n.emplace(std::make_pair(deg, std::move(rec)));

                if (deg == 0) {
                  // this saves some memory since we only need one numerical value
                  // for the constant coefficient
                  saved_num_num[curr_zi_order][ {deg, zi}] = std::make_pair(el.second, sub_count_num);
                }
              }

              if ((int) deg <= curr_deg_num && deg > 0 && curr_zi_order[zi - 2] < deg + 3) {
                saved_num_num[curr_zi_order][ {deg, zi}] = std::make_pair(el.second, sub_count_num);
              }
            }

            for (const auto & el : canonical.second) {
              uint32_t deg = el.first[0];

              if (first_run) {
                if (normalize_to_den && deg == 0)
                  continue;

                PolyReconst rec(n - 1, deg, true);
                coef_d.emplace(std::make_pair(deg, std::move(rec)));

                if (deg == 0)
                  saved_num_den[curr_zi_order][ {deg, zi}] = std::make_pair(el.second, sub_count_den);
              }

              if ((int) deg <= curr_deg_den && deg > 0 && curr_zi_order[zi - 2] < deg + 3) {
                saved_num_den[curr_zi_order][ {deg, zi}] = std::make_pair(el.second, sub_count_den);
              }
            }

            if (first_run) first_run = false;

            uint32_t zi_num = curr_deg_num > 0 ? coef_n[curr_deg_num].get_zi() + 1 : 0;
            uint32_t zi_den = curr_deg_den > 0 ? coef_d[curr_deg_den].get_zi() + 1 : 0;

            // interpolate the numerator
            if (curr_deg_num >= 0) {
              PolyReconst rec = coef_n[curr_deg_num];

              if (rec.get_zi_order() == curr_zi_order) {
                auto res = feed_poly(curr_deg_num, max_deg_num, coef_n, rec,
                                     saved_num_num, sub_num, true);
                curr_deg_num = std::get<0>(res);
                zi_num = std::get<1>(res);
                curr_zi_order_num = std::get<2>(res);
              }
            }

            // interpolate the denominator
            if (curr_deg_den >= 0) {
              PolyReconst rec = coef_d[curr_deg_den];

              // if the numerator is done, get the current zi order of the
              // denominator
              if (curr_deg_num == -1) {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = rec.get_zi_order();
                curr_zi = rec.get_zi() + 1;
                zi = curr_zi;
              }

              if (rec.get_zi_order() == curr_zi_order) {
                auto res = feed_poly(curr_deg_den, max_deg_den, coef_d, rec,
                                     saved_num_den, sub_den, false);
                curr_deg_den = std::get<0>(res);
                zi_den = std::get<1>(res);
                curr_zi_order_den = std::get<2>(res);

                // if the denominator is done, check if the numerator is still undone
                if (curr_deg_den == -1 && curr_deg_num >= 0) {
                  PolyReconst rec_new = coef_n[curr_deg_num];
                  {
                    std::unique_lock<std::mutex> lock(mutex_status);
                    curr_zi_order = rec_new.get_zi_order();
                    curr_zi = rec_new.get_zi() + 1;
                    zi = curr_zi;
                  }

                  if (rec_new.get_zi_order() == std::vector<uint32_t>(curr_zi_order.begin(), curr_zi_order.end())) {
                    auto res_new = feed_poly(curr_deg_num, max_deg_num, coef_n, rec_new,
                                             saved_num_num, sub_num, true);
                    curr_deg_num = std::get<0>(res_new);
                    zi_num = std::get<1>(res_new);
                    curr_zi_order_num = std::get<2>(res_new);
                  }
                }
              }
            }

            // check which poly reconst should be feeded next
            // it is promising that feeding a PolyReconst with a higher
            // zi degree will be finished next leading to less numerical runs
            if (curr_deg_den >= 0 && curr_deg_num >= 0) {
              if (a_grt_b(curr_zi_order_num, curr_zi_order_den)) {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = curr_zi_order_num;
                curr_zi = zi_num;
                zi = zi_num;
              } else {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = curr_zi_order_den;
                curr_zi = zi_den;
                zi = zi_den;
              }
            } else if (curr_deg_num >= 0) {
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = curr_zi_order_num;
              curr_zi = zi_num;
              zi = zi_num;
            } else if (curr_deg_den >= 0) {
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = curr_zi_order_den;
              curr_zi = zi_den;
              zi = zi_den;
            }

            // combine results
            if (curr_deg_den == - 1 && curr_deg_num == -1) {
              saved_num_num.clear();
              saved_num_den.clear();

              first_run = true;

              // Remove normalization due to the shift
              PolynomialFF numerator;
              PolynomialFF denominator;

              for (auto & el : coef_n) {
                PolynomialFF res = el.second.get_result_ff();
                shifted_degs_num.emplace(el.first);

                // TODO: define empty/zero polynomials uniquely
                if (!(res.coefs.size() == 0) && !(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  numerator += res;
                else
                  zero_degs_num.emplace(el.first);
              }

              for (auto & el : coef_d) {
                PolynomialFF res = el.second.get_result_ff();
                shifted_degs_den.emplace(el.first);

                // TODO: define empty/zero polynomials uniquely
                if (!(res.coefs.size() == 0) && !(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  denominator += res;
                else
                  zero_degs_den.emplace(el.first);
              }

              FFInt terminator = 0;
              // check if one can normalize to a single univariate degree, if so
              // set the corresponding degree in equalizer_degree. Check constants
              // first and proceed with denominator/numerator.
              // if the denominator is just a constant, there is no corresponding
              // PolyReconst object. Thus we set the minimal degree to a zero tuple

              if (prime_number == 0) {
                if (const_den != 1) {
                  ff_map dummy_map {};
                  terminator = FFInt(1) - const_den;
                  dummy_map.emplace(std::make_pair(std::vector<uint32_t> (n, 0), terminator));

                  if (normalize_to_den) {
                    denominator += PolynomialFF(n, dummy_map);
                    normalizer_den_num = true;
                  } else {
                    numerator += PolynomialFF(n, dummy_map);
                    normalizer_den_num = false;
                  }

                  normalizer_deg = std::vector<uint32_t> (n, 0);
                } else if (normalize_to_den && numerator.coefs.find(std::vector<uint32_t> (n, 0)) != numerator.coefs.end()) {
                  terminator = numerator.coefs[std::vector<uint32_t> (n, 0)];
                  normalizer_den_num = false;
                  normalizer_deg = std::vector<uint32_t> (n, 0);
                } else if (!normalize_to_den && denominator.coefs.find(std::vector<uint32_t> (n, 0)) != denominator.coefs.end()) {
                  terminator = denominator.coefs[std::vector<uint32_t> (n, 0)];
                  normalizer_den_num = true;
                  normalizer_deg = std::vector<uint32_t> (n, 0);
                } else {
                  for (const auto & el : denominator.coefs) {
                    add_non_solved_den(el.first);

                    if (normalizer_deg.empty())
                      normalizer_deg = el.first;
                    else if (a_grt_b(normalizer_deg, el.first))
                      normalizer_deg = el.first;
                  }

                  for (const auto & candidate : non_solved_degs_den) {
                    if (candidate.second.size() == 1) {
                      terminator = denominator.coefs[candidate.second[0]];
                      normalizer_den_num = true;
                      normalizer_deg = candidate.second[0];
                      break;
                    }
                  }

                  if (terminator == 0) {
                    for (const auto & el : numerator.coefs) {
                      add_non_solved_num(el.first);
                    }

                    for (const auto & candidate : non_solved_degs_num) {
                      if (candidate.second.size() == 1) {
                        terminator = numerator.coefs[candidate.second[0]];
                        normalizer_den_num = false;
                        normalizer_deg = candidate.second[0];
                        break;
                      }
                    }
                  }

                  if (terminator == 0) {
                    // normalize to the minimal degree of the denominator
                    terminator = denominator.coefs[normalizer_deg];
                    normalizer_den_num = true;
                    is_singular_system = true;
                  }
                }

                shifted_max_num_eqn = coef_n.size() + coef_d.size();
              } else {
                if (const_den != 1) {
                  ff_map dummy_map {};
                  terminator = FFInt(1) - const_den;
                  dummy_map.emplace(std::make_pair(std::vector<uint32_t> (n, 0), terminator));

                  if (normalize_to_den) {
                    denominator += PolynomialFF(n, dummy_map);
                    normalizer_den_num = true;
                  } else {
                    numerator += PolynomialFF(n, dummy_map);
                    normalizer_den_num = false;
                  }
                }

                if (normalizer_den_num)
                  terminator = denominator.coefs[normalizer_deg];
                else
                  terminator = numerator.coefs[normalizer_deg];
              }

              curr_zi_order_num = std::vector<uint32_t> (n - 1, 1);
              curr_zi_order_den = std::vector<uint32_t> (n - 1, 1);

              coef_n.clear();
              coef_d.clear();

              // normalize
              FFInt equalizer = FFInt(1) / terminator;

              if (equalizer == 0)
                div_by_zero = true;

              numerator *= equalizer;
              denominator *= equalizer;

              combine_primes(numerator.coefs, denominator.coefs);

              std::unique_lock<std::mutex> lock(mutex_status);
              std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
              curr_zi = 2;
              zi = 2;
            } else if (zi > 0) {
              // set new random
              for (uint32_t tmp_zi = 2; tmp_zi <= n; ++tmp_zi) {
                auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 2]);
                std::unique_lock<std::mutex> lock_statics(mutex_statics);

                if (rand_zi.find(key) == rand_zi.end())
                  rand_zi.emplace(std::make_pair(key, rand_zi[std::make_pair(tmp_zi, 1)].pow(key.second)));
              }
            }

            return;
          } else {
            // build multivariate Vandermonde systems and evaluate them if possible
            if (!is_singular_system) {
              for (const auto & sol : canonical.first) {
                uint32_t key = sol.first[0];

                if (coef_mat_num[key].size() < non_solved_degs_num[key].size())
                  coef_mat_num[key].emplace_back(std::make_pair(sol.second, 0));

                // Solve multivariate Vandermonde system for corresponding degree,
                // remove entry from non_solved_degs and add it to solve_degs
                if (coef_mat_num[key].size() == non_solved_degs_num[key].size()) {
                  solved_num += solve_vandermonde_system(non_solved_degs_num[key], coef_mat_num[key], get_anchor_points());
                  non_solved_degs_num.erase(key);
                  coef_mat_num.erase(key);
                }
              }

              for (const auto & sol : canonical.second) {
                uint32_t key = sol.first[0];

                if (coef_mat_den[key].size() < non_solved_degs_den[key].size())
                  coef_mat_den[key].emplace_back(std::make_pair(sol.second, 0));

                if (coef_mat_den[key].size() == non_solved_degs_den[key].size()) {
                  solved_den += solve_vandermonde_system(non_solved_degs_den[key], coef_mat_den[key], get_anchor_points());
                  non_solved_degs_den.erase(key);
                  coef_mat_den.erase(key);
                }
              }
            } else {
              uint32_t tmp_deg = curr_deg_num;
              int thresh_num = normalize_to_den ? -1 : 0;
              int thresh_den = normalize_to_den ? 0 : -1;

              if (curr_deg_num > thresh_num) {
                uint32_t tmp_deg = curr_deg_num;

                for (uint32_t i = start_deg_num; i <= tmp_deg; ++i) {
                  if (coef_mat_num.find(i) != coef_mat_num.end()) {
                    if (coef_mat_num[i].size() < non_solved_degs_num[i].size())
                      coef_mat_num[i].emplace_back(std::make_pair(canonical.first[ {i}], sub_count_num));

                    if (i == (uint32_t)curr_deg_num && coef_mat_num[i].size() == non_solved_degs_num[i].size()) {
                      set_new_curr_deg_num_singular(i);
                      bool cannot_solve = false;

                      while (!cannot_solve) {
                        std::vector<uint32_t> deleted_degs {};

                        for (const auto & mat : coef_mat_num) {
                          uint32_t tmp_key = mat.first;

                          if (tmp_key == (uint32_t)curr_deg_num && mat.second.size() == non_solved_degs_num[tmp_key].size()) {
                            set_new_curr_deg_num_singular(tmp_key);
                            deleted_degs.emplace_back(tmp_key);
                          }
                        }

                        for (const auto & el : deleted_degs) {
                          coef_mat_num.erase(el);
                        }

                        if (deleted_degs.empty()) cannot_solve = true;

                        coef_mat_num.erase(i);
                      }
                    }
                  }
                }
              }

              if (curr_deg_den > thresh_den) {
                tmp_deg = curr_deg_den;

                for (uint32_t i = start_deg_den; i <= tmp_deg; ++i) {
                  if (coef_mat_den.find(i) != coef_mat_den.end()) {
                    if (coef_mat_den[i].size() < non_solved_degs_den[i].size())
                      coef_mat_den[i].emplace_back(std::make_pair(canonical.second[ {i}], sub_count_den));

                    if (i == (uint32_t)curr_deg_den && coef_mat_den[i].size() == non_solved_degs_den[i].size()) {
                      set_new_curr_deg_den_singular(i);
                      bool cannot_solve = false;

                      while (!cannot_solve) {
                        std::vector<uint32_t> deleted_degs {};

                        for (const auto & mat : coef_mat_den) {
                          uint32_t tmp_key = mat.first;

                          if (tmp_key == (uint32_t)curr_deg_den && mat.second.size() == non_solved_degs_den[tmp_key].size()) {
                            set_new_curr_deg_den_singular(tmp_key);
                            deleted_degs.emplace_back(tmp_key);
                          }
                        }

                        for (const auto & el : deleted_degs) {
                          coef_mat_den.erase(el);
                        }

                        if (deleted_degs.empty()) cannot_solve = true;

                        coef_mat_den.erase(i);
                      }
                    }
                  }
                }
              }
            }

            // promote to next prime and combine results
            if (coef_mat_num.empty() && coef_mat_den.empty()) {
              if (is_singular_system) {
                for (const auto & el : solved_degs_num) solved_num += el.second;

                for (const auto & el : solved_degs_den) solved_den += el.second;

                // normalize
                FFInt terminator = 0;

                if (normalizer_den_num == normalize_to_den &&  normalizer_deg == std::vector<uint32_t> (n)) {
                  if (const_den != 1)
                    terminator = 1 - const_den;
                  else
                    terminator = 1;
                } else {
                  if (normalizer_den_num) {
                    terminator = solved_den.coefs[normalizer_deg];
                  } else {
                    terminator = solved_num.coefs[normalizer_deg];
                  }
                }

                FFInt equalizer = FFInt(1) / terminator;

                for (const auto & el : g_ni) {
                  solved_num.coefs.erase(el.first);
                }

                for (const auto & el : g_di) {
                  solved_den.coefs.erase(el.first);
                }

                solved_num = solved_num * equalizer;
                solved_den = solved_den * equalizer;
              }

              // remove the constant if it is zero
              if (solved_num.coefs.find(std::vector<uint32_t> (n, 0)) != solved_num.coefs.end() && solved_num.coefs[std::vector<uint32_t> (n, 0)] == 0) {
                solved_num.coefs.erase(std::vector<uint32_t> (n, 0));
              }

              if (solved_den.coefs.find(std::vector<uint32_t> (n, 0)) != solved_den.coefs.end() && solved_den.coefs[std::vector<uint32_t> (n, 0)] == 0) {
                solved_den.coefs.erase(std::vector<uint32_t> (n, 0));
              }

              combine_primes(solved_num.coefs, solved_den.coefs);
              {
                std::unique_lock<std::mutex> lock(mutex_status);
                std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
              }
              // reset solved coefficients
              ff_map zero_deg {};
              zero_deg.emplace(std::make_pair(std::vector<uint32_t> (n), 0));
              solved_degs_num = polff_map();
              solved_degs_den = polff_map();
              solved_num.coefs = zero_deg;
              solved_den.coefs = zero_deg;
            } else {
              // increase zi order by 1
              {
                std::unique_lock<std::mutex> lock(mutex_status);
                saved_ti.erase(curr_zi_order);
                std::transform(curr_zi_order.begin(), curr_zi_order.end(),
                curr_zi_order.begin(), [](uint32_t x) {return x + 1;});

                for (uint32_t tmp_zi = 2; tmp_zi <= n; ++tmp_zi) {
                  auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 2]);
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);

                  if (rand_zi.find(key) == rand_zi.end())
                    rand_zi.emplace(std::make_pair(key, rand_zi[std::make_pair(tmp_zi, 1)].pow(key.second)));
                }

                if (!is_singular_system)
                  num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
                else {
                  num_eqn = shifted_max_num_eqn - tmp_solved_coefs_den - tmp_solved_coefs_num;
                }
              }
            }
          }
        }
      } else {
        if (saved_ti.find(feed_zi_ord) == saved_ti.end()) {
          std::vector<std::pair<FFInt, FFInt>> tmp_ti = {std::make_pair(new_ti, num)};
          saved_ti[feed_zi_ord] = tmp_ti;
        } else {
          saved_ti[feed_zi_ord].emplace_back(std::make_pair(new_ti, num));
        }
      }
    }
  }

  std::tuple<int, uint32_t, std::vector<uint32_t>> RatReconst::feed_poly(int curr_deg,
                                                                         uint32_t max_deg, std::unordered_map<uint32_t, PolyReconst>& coef,
  PolyReconst& rec, ff_map_map& saved_num, polff_vec_map& sub_save, bool is_num) {
    uint32_t tmp_zi = rec.get_zi() + 1;
    std::vector<uint32_t> tmp_zi_ord = curr_zi_order;

    while (!rec.is_new_prime()) {
      try {
        std::vector<uint32_t> key = {(uint32_t) curr_deg, tmp_zi};
        auto food_pair = saved_num.at(tmp_zi_ord).at(key);
        FFInt food = food_pair.first;
        uint32_t sub_count = food_pair.second;
        // delete unused saved data
        saved_num[tmp_zi_ord].erase(key);
        // set random values for the yis
        std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);

        // feed to PolyReconst
        // since the constant is just a constant, we do not have to get mutliple
        // numerical values to reconstruct the coefficient
        if (curr_deg == 0) {
          tmp_sol_const_num = 0;
          tmp_sol_const_den = 0;

          while (!rec.is_new_prime()) {
            FFInt sub = 0;

            if (curr_deg != (int)max_deg && sub_count < sub_save[curr_deg].size()) {
              sub = sub_save[curr_deg][sub_count].calc_n_m_1(yis);
            }

            rec.feed(yis, food - sub);
          }
        } else {
          FFInt sub = 0;

          if (curr_deg != (int)max_deg && sub_count < sub_save[curr_deg].size()) {
            sub = sub_save[curr_deg][sub_count].calc_n_m_1(yis);
          }

          rec.feed(yis, food - sub);
          tmp_zi = rec.get_zi() + 1;
          tmp_zi_ord = rec.get_zi_order();
        }
      } catch (std::out_of_range& e) {
        coef[curr_deg] = rec;
        return std::make_tuple(curr_deg, tmp_zi, tmp_zi_ord);
      }

      /*
       * Check for new prime & save subtraction term
       */
      if (rec.is_new_prime()) {
        coef[curr_deg] = rec;

        if (curr_deg > 0) {
          PolynomialFF res = rec.get_result_ff();

          if (!res.zero()) {
            std::vector<uint32_t> zero_deg(n);
            PolynomialFF zero_poly(n, {{zero_deg, 0}});

            //todo only save needed shifts
            for (int i = 0; i < curr_deg; ++i) {
              if (sub_save[(uint32_t)i].size() == 0)
                sub_save[(uint32_t)i] = {zero_poly};
              else
                sub_save[(uint32_t)i].emplace_back(zero_poly);
            }

            PolynomialFF res = rec.get_result_ff();

            std::vector<FFInt> tmp_shift;
            {
              std::unique_lock<std::mutex> lock_statics(mutex_statics);
              tmp_shift = shift;
            }

            PolynomialFF sub_pol = rec.get_result_ff().add_shift(tmp_shift);

            for (auto & el : sub_pol.coefs) {
              int tmp_deg = 0;

              for (const auto & deg : el.first) tmp_deg += deg;

              if (tmp_deg < curr_deg) {
                for (auto & tmp_sub : sub_save[(uint32_t)tmp_deg]) {
                  tmp_sub += PolynomialFF(n, {{el.first, el.second}});
                }
              }
            }

            if (!is_num) {
              sub_count_den ++;

              if (normalize_to_den) {
                std::vector<FFInt> tmp_yis(n - 1, 0);
                const_den += sub_save[0].back().calc_n_m_1(tmp_yis);
              }
            } else {
              sub_count_num ++;

              if (!normalize_to_den) {
                std::vector<FFInt> tmp_yis(n - 1, 0);
                const_den += sub_save[0].back().calc_n_m_1(tmp_yis);
              }
            }
          }
        }

        sub_save[curr_deg] = std::vector<PolynomialFF>();

        /*
         * Remove already solved coefficients from Gauss eliminiation
         */
        curr_deg--;
        bool found = false;

        if (is_num) {
          if (normalize_to_den) {
            if (curr_deg > -1) {
              while (!found) {
                if (coef_n.find(curr_deg) == coef_n.end())
                  curr_deg --;
                else
                  found = true;

                if (curr_deg == -1)
                  found = true;
              }
            }
          } else {
            if (curr_deg > 0) {
              while (!found) {
                if (coef_n.find(curr_deg) == coef_n.end())
                  curr_deg --;
                else
                  found = true;

                if (curr_deg == 0) {
                  curr_deg = -1;
                  found = true;
                }
              }
            }

            if (curr_deg == 0)
              curr_deg = -1;
          }

          tmp_solved_coefs_num ++;
        } else {
          if (normalize_to_den) {
            if (curr_deg > 0) {
              while (!found) {
                if (coef_d.find(curr_deg) == coef_d.end())
                  curr_deg --;
                else
                  found = true;

                if (curr_deg == 0) {
                  curr_deg = -1;
                  found = true;
                }
              }
            }

            if (curr_deg == 0)
              curr_deg = -1;
          } else {
            if (curr_deg > -1) {
              while (!found) {
                if (coef_d.find(curr_deg) == coef_d.end())
                  curr_deg --;
                else
                  found = true;

                if (curr_deg == -1)
                  found = true;
              }
            }
          }

          tmp_solved_coefs_den ++;
        }

        {
          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = max_deg_den + max_deg_num + 1
                    - tmp_solved_coefs_num - tmp_solved_coefs_den
                    - tmp_sol_const_num - tmp_sol_const_den;
        }

        if (curr_deg >= 0) {
          rec = coef[curr_deg];
          std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end(), 1);
          tmp_zi = rec.get_zi() + 1;
        } else break;
      }
    }

    tmp_zi = rec.get_zi() + 1;
    return std::make_tuple(curr_deg, tmp_zi, tmp_zi_ord);
  }

  void RatReconst::combine_primes(ff_map& numerator, ff_map& denominator) {
    std::vector<uint32_t> tmp_deg_num {};
    std::vector<uint32_t> tmp_deg_den {};

    sub_count_den = 0;
    sub_count_num = 0;
    saved_ti = ff_vec_map();
    tmp_sol_const_num = 0;
    tmp_sol_const_den = 0;

    if (is_singular_system && !div_by_zero) {
      for (const auto & el : non_solved_degs_num) {
        tmp_deg_num.emplace_back(el.first);
      }

      for (const auto & el : non_solved_degs_den) {
        tmp_deg_den.emplace_back(el.first);
      }
    }

    non_solved_degs_den.clear();
    non_solved_degs_num.clear();
    sub_num = polff_vec_map();
    sub_den = polff_vec_map();
    saved_num_num = ff_map_map();
    saved_num_den = ff_map_map();
    coef_d = std::unordered_map<uint32_t, PolyReconst>();
    coef_n = std::unordered_map<uint32_t, PolyReconst>();

    if (!div_by_zero) {
      if (!use_chinese_remainder) {
        combined_ni = convert_to_mpz(numerator);
        combined_di = convert_to_mpz(denominator);

        mpz_map combined_ni_back = combined_ni;

        for (auto & c_ni : combined_ni_back) {

          if (!normalizer_den_num && c_ni.first == normalizer_deg)
            remove_ni(c_ni.first, RationalNumber(1, 1));
          else
            add_non_solved_num(c_ni.first);
        }

        mpz_map combined_di_back = combined_di;

        for (auto & c_di : combined_di_back) {

          if (normalizer_den_num && c_di.first == normalizer_deg)
            remove_di(c_di.first, RationalNumber(1, 1));
          else
            add_non_solved_den(c_di.first);
        }


        if (is_singular_system) {
          check_for_solved_degs(tmp_deg_num, true);

          if (is_singular_system)
            check_for_solved_degs(tmp_deg_den, false);
        }
      } else {
        for (const auto & el : numerator) {
          if (combined_ni.find(el.first) == combined_ni.end() && g_ni.find(el.first) == g_ni.end())
            combined_ni.emplace(std::make_pair(el.first, 0));
        }

        for (const auto & el : denominator) {
          if (combined_di.find(el.first) == combined_di.end() && g_di.find(el.first) == g_di.end())
            combined_di.emplace(std::make_pair(el.first, 0));
        }

        for (const auto & el : combined_ni) {
          if (numerator.find(el.first) == numerator.end() && g_ni.find(el.first) == g_ni.end())
            numerator.emplace(std::make_pair(el.first, 0));
        }

        for (const auto & el : combined_di) {
          if (denominator.find(el.first) == denominator.end() && g_di.find(el.first) == g_di.end())
            denominator.emplace(std::make_pair(el.first, 0));
        }

        mpz_map combined_ni_back = combined_ni;
        mpz_map combined_di_back = combined_di;
        mpz_class combined_prime_back = combined_prime;
        std::pair<mpz_class, mpz_class> p1;
        std::pair<mpz_class, mpz_class> p2;
        std::pair<mpz_class, mpz_class> p3;

        //numerator
        for (auto it = combined_ni.begin(); it != combined_ni.end(); ++it) {
          p1 = std::make_pair(it->second, combined_prime);
          p2 = std::make_pair(mpz_class(numerator[it->first].n), FFInt::p);
          p3 = run_chinese_remainder(p1, p2);
          combined_ni[it->first] = p3.first;
        }

        // denominator
        for (auto it = combined_di.begin(); it != combined_di.end(); ++it) {
          p1 = std::make_pair(it->second, combined_prime);
          p2 = std::make_pair(mpz_class(denominator[it->first].n), FFInt::p);
          p3 = run_chinese_remainder(p1, p2);
          combined_di[it->first] = p3.first;
        }

        combined_prime = p3.second;

        // Remove already known coefficients from solve algorithm to save numerical runs
        for (auto & c_ni : combined_ni_back) {
          if (g_ni.find(c_ni.first) == g_ni.end()) {
            RationalNumber last_rn_wang;
            RationalNumber curr_rn_wang;
            bool last_wang;
            bool curr_wang;
            auto res = get_rational_coef(c_ni.second, combined_prime_back);
            last_wang = res.first;

            if (last_wang)
              last_rn_wang = res.second;

            res = get_rational_coef(combined_ni[c_ni.first], combined_prime);
            curr_wang = res.first;

            if (curr_wang)
              curr_rn_wang = res.second;


            RationalNumber last_rn_monagan;
            RationalNumber curr_rn_monagan;
            bool last_monagan;
            bool curr_monagan;

            res = get_rational_coef_mqrr(c_ni.second, combined_prime_back);
            last_monagan = res.first;

            if (last_monagan)
              last_rn_monagan = res.second;

            res = get_rational_coef_mqrr(combined_ni[c_ni.first], combined_prime);
            curr_monagan = res.first;

            if (curr_monagan)
              curr_rn_monagan = res.second;


            if (last_wang && curr_wang && last_rn_wang == curr_rn_wang) {
              remove_ni(c_ni.first, curr_rn_wang);
              continue;
            }

            if (last_monagan && curr_monagan && last_rn_monagan == curr_rn_monagan) {
              remove_ni(c_ni.first, curr_rn_monagan);
              continue;
            }

            if (c_ni.second == combined_ni[c_ni.first]) {
              RationalNumber rn = RationalNumber(c_ni.second, 1);
              remove_ni(c_ni.first, rn);
            } else
              add_non_solved_num(c_ni.first);
          }
        }

        for (auto & c_di : combined_di_back) {
          if (g_di.find(c_di.first) == g_di.end()) {
            RationalNumber last_rn_wang;
            RationalNumber curr_rn_wang;
            bool last_wang;
            bool curr_wang;
            auto res = get_rational_coef(c_di.second, combined_prime_back);
            last_wang = res.first;

            if (last_wang)
              last_rn_wang = res.second;

            res = get_rational_coef(combined_di[c_di.first], combined_prime);
            curr_wang = res.first;

            if (curr_wang)
              curr_rn_wang = res.second;


            RationalNumber last_rn_monagan;
            RationalNumber curr_rn_monagan;
            bool last_monagan;
            bool curr_monagan;

            res = get_rational_coef_mqrr(c_di.second, combined_prime_back);
            last_monagan = res.first;

            if (last_monagan)
              last_rn_monagan = res.second;

            res = get_rational_coef_mqrr(combined_di[c_di.first], combined_prime);
            curr_monagan = res.first;

            if (curr_monagan)
              curr_rn_monagan = res.second;

            if (last_wang && curr_wang && last_rn_wang == curr_rn_wang) {
              remove_di(c_di.first, curr_rn_wang);
              continue;
            }

            if (last_monagan && curr_monagan && last_rn_monagan == curr_rn_monagan) {
              remove_di(c_di.first, curr_rn_monagan);
              continue;
            }

            if (c_di.second == combined_di[c_di.first]) {
              RationalNumber rn = RationalNumber(c_di.second, 1);
              remove_di(c_di.first, rn);
            } else
              add_non_solved_den(c_di.first);
          }
        }

        if (is_singular_system) {
          check_for_solved_degs(tmp_deg_num, true);

          if (is_singular_system)
            check_for_solved_degs(tmp_deg_den, false);
        }
      }

      if (prime_number + 1 >= interpolations) {
        if (is_singular_system) {
          {
            std::unique_lock<std::mutex> lock_statics(mutex_statics);
            need_prime_shift = true;
          }

          for (const auto & el : g_ni) {
            if (normalize_to_den)
              add_non_solved_num(el.first);
            else if (el.first != std::vector<uint32_t> (n))
              add_non_solved_num(el.first);
          }

          for (const auto & el : g_di) {
            if (!normalize_to_den)
              add_non_solved_den(el.first);
            else if (el.first != std::vector<uint32_t> (n))
              add_non_solved_den(el.first);
          }

          if (normalize_to_den) {
            curr_deg_num = max_deg_num;

            if (max_deg_den > 0)
              curr_deg_den = max_deg_den;
          } else {
            curr_deg_den = max_deg_den;

            if (max_deg_num > 0)
              curr_deg_num = max_deg_num;
          }

          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = shifted_max_num_eqn;
        } else {
          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
        }
      }

      for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

      for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};
    } else {
      div_by_zero = false;
      interpolations ++;
    }

    if (prime_number + 1 < interpolations) {
      max_deg_num = -1;
      first_run = true;
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      need_prime_shift = true;
    }

    const_den = 0;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      ++prime_number;
      queue = std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>>();
      saved_ti.clear();
      new_prime = true;
    }
    tmp_solved_coefs_den = 0;
    tmp_solved_coefs_num = 0;

    // Check if the state should be written out after this prime
    if (tag.size() > 0)
      save_state();
  }

  RationalFunction RatReconst::get_result() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (done) {
      if (result.numerator.coefs.empty()) {
        Polynomial numerator;
        Polynomial denominator;

        numerator = Polynomial(g_ni);
        denominator = Polynomial(g_di);
        g_ni.clear();
        g_di.clear();
        //numerator.sort();
        //denominator.sort();

        result = RationalFunction(numerator, denominator);
      }

      return result;
    } else {
      ERROR_MSG("Trying to access unfinished result.");
      std::exit(-1);
    }
  }

  bool RatReconst::rec_rat_coef() {
    bool run_test = true;
    std::vector<const std::vector<uint32_t>*> promoted_n;
    std::vector<const std::vector<uint32_t>*> promoted_d;

    for (const auto & ci : combined_ni) {
      mpz_class a = ci.second;
      auto res = get_rational_coef(a, combined_prime);

      if (res.first) {
        g_ni[ci.first] = res.second;
        promoted_n.emplace_back(&ci.first);
      } else {
        res = get_rational_coef_mqrr(a, combined_prime);

        if (res.first) {
          g_ni[ci.first] = res.second;
          promoted_n.emplace_back(&ci.first);
        } else {
          run_test = false;
          break;
        }
      }
    }

    if (run_test) {
      for (const auto & ci : combined_di) {
        mpz_class a = ci.second;
        auto res = get_rational_coef(a, combined_prime);

        if (res.first) {
          g_di[ci.first] = res.second;
          promoted_d.emplace_back(&ci.first);
        } else {
          res = get_rational_coef_mqrr(a, combined_prime);

          if (res.first) {
            g_di[ci.first] = res.second;
            promoted_d.emplace_back(&ci.first);
          } else {
            run_test = false;
            break;
          }
        }
      }
    }

    if (!run_test) {
      for (const auto & ci : promoted_n) g_ni.erase(*ci);

      for (const auto & ci : promoted_d) g_di.erase(*ci);
    }

    return run_test;
  }

  std::pair<ff_map, ff_map> RatReconst::solve_gauss() {
    std::vector<FFInt> results = solve_gauss_system(coef_mat, num_eqn);
    coef_mat.clear();
    num_sub_den.clear();
    num_sub_num.clear();

    ff_map numerator;
    ff_map denominator;
    uint32_t counter = 0;

    for (const auto & el : coef_n) {
      uint32_t tmp_key = el.first;

      if (tmp_key == 0 && tmp_sol_const_num == 1)
        continue;

      if ((int)tmp_key <= curr_deg_num) {
        std::vector<uint32_t> power = {tmp_key};
        numerator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }
    }

    for (const auto & el : coef_d) {
      uint32_t tmp_key = el.first;

      if (tmp_key == 0 && tmp_sol_const_den == 1)
        continue;

      if ((int) tmp_key <= curr_deg_den) {
        std::vector<uint32_t> power = {tmp_key};
        denominator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }
    }

    return std::make_pair(numerator, denominator);
  }

  std::pair<ff_map, ff_map> RatReconst::solve_homogenized_multi_gauss() {
    std::vector<FFInt> results = solve_gauss_system(coef_mat, num_eqn);
    coef_mat.clear();
    num_sub_den.clear();
    num_sub_num.clear();

    ff_map numerator;
    ff_map denominator;
    uint32_t counter = 0;

    if (!is_singular_system) {
      for (const auto & el : non_solved_degs_num) {
        std::vector<uint32_t> power = {el.first};
        numerator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }

      for (const auto & el : non_solved_degs_den) {
        std::vector<uint32_t> power = {el.first};
        denominator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }
    } else {
      for (const auto & el : shifted_degs_num) {
        uint32_t tmp_key = el;

        if ((int)tmp_key <= curr_deg_num) {
          std::vector<uint32_t> power = {tmp_key};
          numerator.emplace(std::make_pair(std::move(power), results[counter]));
          counter ++;
        }
      }

      for (const auto & el : shifted_degs_den) {
        uint32_t tmp_key = el;

        if ((int) tmp_key <= curr_deg_den) {
          std::vector<uint32_t> power = {tmp_key};
          denominator.emplace(std::make_pair(std::move(power), results[counter]));
          counter ++;
        }
      }
    }

    return std::make_pair(numerator, denominator);
  }

  RationalFunction RatReconst::normalize(RationalFunction& rf) {
    RationalNumber equalizer = rf.denominator.coefs[0].coef;
    RationalNumber terminator(equalizer.denominator, equalizer.numerator);

    rf.numerator *= terminator;
    rf.denominator *= terminator;
    return rf;
  }

  bool RatReconst::test_guess(const FFInt& num, const FFInt& ti) {
    ff_map g_ff_ni = convert_to_ffint(g_ni);
    ff_map g_ff_di = convert_to_ffint(g_di);
    PolynomialFF g_ny(n, g_ff_ni);
    PolynomialFF g_dy(n, g_ff_di);
    std::vector<FFInt> yis = std::vector<FFInt> (n);
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yis[0] = ti + shift[0];

      for (uint32_t i = 1; i < n; ++i) {
        yis[i] = ti * rand_zi[std::make_pair(i + 1, 1)].pow(curr_zi_order[i - 1]) + shift[i];
      }
    }

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;
  }

  void RatReconst::remove_ni(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn) {
    g_ni[deg_vec] =  rn;
    combined_ni.erase(deg_vec);
  }

  void RatReconst::remove_di(const std::vector<uint32_t>& deg_vec, const RationalNumber& rn) {
    g_di[deg_vec] =  rn;
    combined_di.erase(deg_vec);
  }

  void RatReconst::disable_shift() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    shift = std::vector<FFInt> (n, 0);
  }

  void RatReconst::build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);

    for (uint32_t i = 0; i < n; ++i) {
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        yis[i] = yis[i] * tmp_ti + shift[i];
      }
    }

    FFInt res;

    if (normalize_to_den)
      res = (1 - const_den) * tmp_num;
    else
      res = -(1 - const_den);

    for (auto & el : coef_n) {
      uint32_t tmp_key = el.first;

      if (tmp_key == 0 && tmp_sol_const_num == 1) {
        res -= saved_num_num.at(std::vector<uint32_t> (n - 1, 1)).at( {0, 2}).first;

        if (curr_deg_num < max_deg_num) {
          res += sub_num[0].front().calc(yis);
        }
      } else if ((int) tmp_key > curr_deg_num)
        res -= el.second.get_result_ff().calc(yis);
      else
        eq.emplace_back(tmp_ti.pow(tmp_key));
    }

    for (auto & el : coef_d) {
      uint32_t tmp_key = el.first;

      if (tmp_key == 0 && tmp_sol_const_den == 1) {
        res += saved_num_den.at(std::vector<uint32_t> (n - 1, 1)).at( {0, 2}).first* tmp_num;

        if (curr_deg_den < max_deg_den)
          res -= sub_den[0].front().calc(yis) * tmp_num;
      } else if ((int) tmp_key > curr_deg_den)
        res += el.second.get_result_ff().calc(yis) * tmp_num;
      else
        eq.emplace_back(-tmp_ti.pow(tmp_key) * tmp_num);
    }

    eq.emplace_back(res);

    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);

    if (!is_singular_system) {
      // Build system of equations; in combined_.. are the non-solved coefficients
      for (const auto & pow_vec : non_solved_degs_num) {
        eq.emplace_back(tmp_ti.pow(pow_vec.first));
      }

      for (const auto & pow_vec : non_solved_degs_den) {
        eq.emplace_back(FFInt(0) - tmp_num * tmp_ti.pow(pow_vec.first));
      }

      // Build result vector including subtracted coefficients which have already
      // been solved
      eq.emplace_back(0);

      if (coef_mat.size() == 0) {
        yis.erase(yis.begin());
        num_sub_num = solved_num.calc_n_m_1_map(yis);
        num_sub_den = solved_den.calc_n_m_1_map(yis);
      }

      for (const auto & el : num_sub_num) {
        eq.back() -= el.second * tmp_ti.pow(el.first);
      }

      for (const auto & el : num_sub_den) {
        eq.back() += el.second * tmp_ti.pow(el.first) * tmp_num;
      }

      coef_mat.emplace_back(std::move(eq));
    } else {
      FFInt res;

      for (uint32_t i = 0; i < n; ++i) {
        {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          yis[i] = yis[i] * tmp_ti + shift[i];
        }
      }

      if (normalize_to_den)
        res = (1 - const_den) * tmp_num;
      else
        res = -(1 - const_den);

      for (const auto & el : shifted_degs_num) {
        uint32_t tmp_key = el;

        if ((int) tmp_key > curr_deg_num)
          res -= solved_degs_num[tmp_key].calc(yis);
        else
          eq.emplace_back(tmp_ti.pow(tmp_key));
      }

      for (const auto & el : shifted_degs_den) {
        uint32_t tmp_key = el;

        if ((int) tmp_key > curr_deg_den)
          res += solved_degs_den[tmp_key].calc(yis) * tmp_num;
        else
          eq.emplace_back(-tmp_ti.pow(tmp_key) * tmp_num);
      }

      eq.emplace_back(res);

      coef_mat.emplace_back(std::move(eq));
    }
  }

  void RatReconst::generate_anchor_points() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    rand_zi.clear();

    for (uint32_t tmp_zi = 2; tmp_zi <= n; ++tmp_zi) {
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 0), 1));
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 1), FFInt(xorshift64star())));
    }

    PolyReconst rec(n - 1, 0, true);
    std::vector<FFInt> anchor_points {};

    for (uint32_t tmp_zi = 2; tmp_zi <= n; ++tmp_zi) {
      anchor_points.emplace_back(rand_zi[std::make_pair(tmp_zi, 1)]);
    }

    rec.set_anchor_points(anchor_points, true);
  }

  void RatReconst::add_non_solved_num(const std::vector<uint32_t>& deg) {
    uint32_t degree = 0;

    for (const auto & el : deg) degree += el;

    non_solved_degs_num[degree].emplace_back(deg);
  }

  void RatReconst::add_non_solved_den(const std::vector<uint32_t>& deg) {
    uint32_t degree = 0;

    for (const auto & el : deg) degree += el;

    non_solved_degs_den[degree].emplace_back(deg);
  }

  void RatReconst::check_for_solved_degs(const std::vector<uint32_t>& uni_degs, const bool is_num) {
    for (const auto & el : uni_degs) {
      if (is_num) {
        if (non_solved_degs_num.find(el) == non_solved_degs_num.end()) {
          is_singular_system = false;
          break;
        }
      } else {
        if (non_solved_degs_den.find(el) == non_solved_degs_den.end()) {
          is_singular_system = false;
          break;
        }
      }
    }
  }

  FFInt RatReconst::get_rand_zi(uint32_t zi, uint32_t order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return rand_zi.at(std::make_pair(zi, order));
  }

  std::vector<FFInt> RatReconst::get_rand_zi_vec(const std::vector<uint32_t>& order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    std::vector<FFInt> res {};

    for (uint32_t i = 2; i <= n; ++i) {
      res.emplace_back(rand_zi.at(std::make_pair(i, order[i - 2])));
    }

    return res;
  }

  FFInt RatReconst::get_zi_shift(uint32_t zi) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return shift[zi - 1];
  }

  bool RatReconst::is_shift_working() {
    std::unique_lock<std::mutex> lock_status(mutex_status);
    std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
    done = false;
    zi = 1;
    return shift_works;
  }

  void RatReconst::accept_shift() {
    scan = false;
  }

  bool RatReconst::get_is_interpolating() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return is_interpolating;
  }

  std::vector<FFInt> RatReconst::get_zi_shift_vec() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return shift;
  }

  bool RatReconst::need_shift() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    bool tmp = need_prime_shift;
    set_singular_system = need_prime_shift;
    need_prime_shift = false;
    return tmp;
  }

  void RatReconst::set_new_curr_deg_num_singular(uint32_t key) {
    if (curr_deg_num < max_deg_num) {
      for (uint32_t i = 0; i < coef_mat_num[key].size(); ++i) {
        auto tmp_pair = coef_mat_num[key][i];

        if (tmp_pair.second < sub_num[key].size()) {
          std::vector<uint32_t> tmp_zi_ord(n - 1, i + 1);
          std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
          tmp_pair.first -= sub_num[key][tmp_pair.second].calc_n_m_1(yis);
          coef_mat_num[key][i] = tmp_pair;
        }
      }
    }

    solved_degs_num[key] = solve_vandermonde_system(non_solved_degs_num[key], coef_mat_num[key], get_anchor_points());

    std::vector<uint32_t> zero_deg(n);
    PolynomialFF zero_poly(n, {{zero_deg, 0}});

    for (const auto & el : coef_mat_num) {
      uint32_t tmp_key = el.first;

      if (sub_num[tmp_key].size() == 0)
        sub_num[tmp_key] = {zero_poly};
      else
        sub_num[tmp_key].emplace_back(zero_poly);
    }

    if (!normalize_to_den) {
      if (sub_num[0].size() == 0)
        sub_num[0] = {zero_poly};
      else
        sub_num[0].emplace_back(zero_poly);
    }

    if (key > 0) {
      std::vector<FFInt> tmp_shift;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        tmp_shift = shift;
      }
      PolynomialFF sub_pol = solved_degs_num[key].add_shift(tmp_shift);

      for (auto & el : sub_pol.coefs) {
        int tmp_deg = 0;

        for (const auto & deg : el.first) tmp_deg += deg;

        if (tmp_deg < curr_deg_num && (coef_mat_num.find(tmp_deg) != coef_mat_num.end() || tmp_deg == 0)) {
          for (auto & tmp_sub : sub_num[(uint32_t)tmp_deg]) {
            tmp_sub += PolynomialFF(n, {{el.first, el.second}});
          }
        }
      }
    }

    sub_num[key] = std::vector<PolynomialFF>();
    sub_count_num ++;
    bool found = false;
    curr_deg_num --;

    if (normalize_to_den) {
      if (curr_deg_num > -1) {
        while (!found) {
          if (coef_mat_num.find(curr_deg_num) == coef_mat_num.end()) {
            if (zero_degs_num.find(curr_deg_num) != zero_degs_num.end())
              tmp_solved_coefs_num ++;

            curr_deg_num --;
          } else
            found = true;

          if (curr_deg_num == -1)
            found = true;
        }
      }
    } else {
      std::vector<FFInt> tmp_yis(n - 1, 0);
      const_den += sub_num[0].back().calc_n_m_1(tmp_yis);

      if (curr_deg_num > 0) {
        while (!found) {
          if (coef_mat_num.find(curr_deg_num) == coef_mat_num.end()) {
            if (zero_degs_num.find(curr_deg_num) != zero_degs_num.end())
              tmp_solved_coefs_num ++;

            curr_deg_num --;
          } else
            found = true;

          if (curr_deg_num == 0) {
            curr_deg_num = -1;
            found = true;
          }
        }
      }

      if (curr_deg_num == 0)
        curr_deg_num = -1;
    }

    tmp_solved_coefs_num ++;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = shifted_max_num_eqn - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }
  }

  void RatReconst::set_new_curr_deg_den_singular(uint32_t key) {
    if (curr_deg_den < max_deg_den) {
      for (uint32_t i = 0; i < coef_mat_den[key].size(); ++i) {
        auto tmp_pair = coef_mat_den[key][i];

        if (tmp_pair.second < sub_den[key].size()) {
          std::vector<uint32_t> tmp_zi_ord(n - 1, i + 1);
          std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
          tmp_pair.first -= sub_den[key][tmp_pair.second].calc_n_m_1(yis);
          coef_mat_den[key][i] = tmp_pair;
        }
      }
    }

    solved_degs_den[key] = solve_vandermonde_system(non_solved_degs_den[key], coef_mat_den[key], get_anchor_points());

    std::vector<uint32_t> zero_deg(n);
    PolynomialFF zero_poly(n, {{zero_deg, 0}});


    for (const auto & el : coef_mat_den) {
      uint32_t tmp_key = el.first;

      if (sub_den[tmp_key].size() == 0)
        sub_den[tmp_key] = {zero_poly};
      else
        sub_den[tmp_key].emplace_back(zero_poly);
    }

    if (normalize_to_den) {
      if (sub_den[0].size() == 0)
        sub_den[0] = {zero_poly};
      else
        sub_den[0].emplace_back(zero_poly);
    }

    if (curr_deg_den > 0) {
      std::vector<FFInt> tmp_shift;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        tmp_shift = shift;
      }
      PolynomialFF sub_pol = solved_degs_den[key].add_shift(tmp_shift);

      for (auto & el : sub_pol.coefs) {
        int tmp_deg = 0;

        for (const auto & deg : el.first) tmp_deg += deg;

        if (tmp_deg < curr_deg_den && (coef_mat_den.find(tmp_deg) != coef_mat_den.end() || tmp_deg == 0)) {
          for (auto & tmp_sub : sub_den[(uint32_t)tmp_deg]) {
            tmp_sub += PolynomialFF(n, {{el.first, el.second}});
          }
        }
      }
    }

    sub_den[key] = std::vector<PolynomialFF>();
    sub_count_den ++;

    bool found = false;
    curr_deg_den --;

    if (normalize_to_den) {
      std::vector<FFInt> tmp_yis(n - 1, 0);
      const_den += sub_den[0].back().calc_n_m_1(tmp_yis);

      if (curr_deg_den > 0) {
        while (!found) {
          if (coef_mat_den.find(curr_deg_den) == coef_mat_den.end()) {
            if (zero_degs_den.find(curr_deg_den) != zero_degs_den.end())
              tmp_solved_coefs_den ++;

            curr_deg_den --;
          } else
            found = true;

          if (curr_deg_den == 0) {
            curr_deg_den = -1;
            found = true;
          }
        }
      }

      if (curr_deg_den == 0)
        curr_deg_den = -1;
    } else {
      if (curr_deg_den > -1) {
        while (!found) {
          if (coef_mat_den.find(curr_deg_den) == coef_mat_den.end()) {
            if (zero_degs_den.find(curr_deg_den) != zero_degs_den.end())
              tmp_solved_coefs_den ++;

            curr_deg_den --;
          } else
            found = true;

          if (curr_deg_den == -1)
            found = true;
        }
      }
    }

    tmp_solved_coefs_den ++;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = shifted_max_num_eqn - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }
  }

  void RatReconst::set_tag(const std::string& tag_) {
    tag = tag_;
  }

  void RatReconst::save_state() {
    mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::ofstream file;
    std::string file_name = "ff_save/" + tag + "_" + std::to_string(prime_number) + ".txt";
    file.open(file_name.c_str());
    file << "combined_prime\n" << combined_prime.get_str() << "\n";
    file << "is_done\n" << is_done() << "\n";
    file << "max_deg_num\n" << max_deg_num << "\n";
    file << "max_deg_den\n" << max_deg_den << "\n";
    file << "need_prime_shift\n" << need_prime_shift << "\n";
    file << "normalizer_deg\n";
    std::string tmp_vec = "";

    for (const auto & deg : normalizer_deg) {
      tmp_vec += std::to_string(deg) + " ";
    }

    tmp_vec.substr(0, tmp_vec.size() - 1);
    tmp_vec += std::string("\n");
    file << tmp_vec;

    file << "normalize_to_den\n" << normalize_to_den << "\n";
    file << "normalizer_den_num\n" << normalizer_den_num << "\n";
    file << "shifted_max_num_eqn\n" << shifted_max_num_eqn << "\n";
    file << "shift\n";

    std::string tmp_shift = "";

    {
      std::unique_lock<std::mutex> lock(mutex_statics);

      for (const auto & sft : shift) {
        if (sft > 0)
          tmp_shift += "1 ";
        else
          tmp_shift += "0 ";
      }

    }
    tmp_shift.substr(0, tmp_shift.size() - 1);
    tmp_shift += std::string("\n");

    file << tmp_shift;

    file << "shifted_degs_num\n";
    std::string tmp_str = "";

    for (const auto & el : shifted_degs_num) {
      tmp_str += std::to_string(el) + " ";
    }

    if (shifted_degs_num.size() > 0) {
      tmp_str.substr(0, tmp_str.size() - 1);
      tmp_str += std::string("\n");
      file << tmp_str;
    }

    file << "shifted_degs_den\n";
    tmp_str = "";

    for (const auto & el : shifted_degs_den) {
      tmp_str += std::to_string(el) + " ";
    }

    if (shifted_degs_den.size() > 0) {
      tmp_str.substr(0, tmp_str.size() - 1);
      tmp_str += std::string("\n");
      file << tmp_str;
    }

    file << "zero_degs_num\n";
    tmp_str = "";

    for (const auto & el : zero_degs_num) {
      tmp_str += std::to_string(el) + " ";
    }

    if (zero_degs_num.size() > 0) {
      tmp_str.substr(0, tmp_str.size() - 1);
      tmp_str += std::string("\n");
      file << tmp_str;
    }

    file << "zero_degs_den\n";
    tmp_str = "";

    for (const auto & el : zero_degs_den) {
      tmp_str += std::to_string(el) + " ";
    }

    if (zero_degs_den.size() > 0) {
      tmp_str.substr(0, tmp_str.size() - 1);
      tmp_str += std::string("\n");
      file << tmp_str;
    }

    file << "g_ni\n";

    for (const auto & el : g_ni) {
      std::string tmp_entry = "";

      for (const auto & deg : el.first) tmp_entry += std::to_string(deg) + " ";

      tmp_entry += std::string(el.second.numerator.get_str()) + " " + std::string(el.second.denominator.get_str()) + std::string("\n");
      file << tmp_entry;
    }

    file << "g_di\n";

    for (const auto & el : g_di) {
      std::string tmp_entry = "";

      for (const auto & deg : el.first) tmp_entry += std::to_string(deg) + std::string(" ");

      tmp_entry += std::string(el.second.numerator.get_str()) + " " + std::string(el.second.denominator.get_str()) + std::string("\n");
      file << tmp_entry;
    }

    file << "combined_ni\n";

    for (const auto & el : combined_ni) {
      std::string tmp_entry = "";

      for (const auto & deg : el.first) tmp_entry += std::to_string(deg) + std::string(" ");

      tmp_entry += std::string(el.second.get_str()) + std::string("\n");
      file << tmp_entry;
    }

    file << "combined_di\n";

    for (const auto & el : combined_di) {
      std::string tmp_entry = "";

      for (const auto & deg : el.first) tmp_entry += std::to_string(deg) + std::string(" ");

      tmp_entry += std::string(el.second.get_str()) + std::string("\n");
      file << tmp_entry;
    }

    file << "interpolations\n" << interpolations << "\n";

    file.close();

    if (prime_number > 0) {
      std::string old_file_name = "ff_save/" + tag + "_" + std::to_string(prime_number - 1) + ".txt";

      if (std::remove(old_file_name.c_str()) != 0)
        WARNING_MSG("The previously saved file could not be deleted.");
    }
  }

  void RatReconst::start_from_saved_file(std::string file_name) {
    std::string line;
    std::ifstream file(file_name.c_str());
    bool first = true;
    parse_prime_number(file_name);

    if (file.is_open()) {
      while (std::getline(file, line)) {
        if (first) {
          first = false;

          if (line != "combined_prime") {
            ERROR_MSG("Wrong input format! Has to start with 'combined_prime'!");
            std::exit(-1);
          }

          curr_parsed_variable = COMBINED_PRIME;
          parsed_variables[COMBINED_PRIME] = true;
        } else {
          if (line == "is_done") {
            curr_parsed_variable = IS_DONE;
            parsed_variables[IS_DONE] = true;
          } else if (line == "max_deg_num") {
            curr_parsed_variable = MAX_DEG_NUM;
            parsed_variables[MAX_DEG_NUM] = true;
          } else if (line == "max_deg_den") {
            curr_parsed_variable = MAX_DEG_DEN;
            parsed_variables[MAX_DEG_DEN] = true;
          } else if (line == "need_prime_shift") {
            curr_parsed_variable = NEED_PRIME_SHIFT;
            parsed_variables[NEED_PRIME_SHIFT] = true;
          } else if (line == "normalizer_deg") {
            curr_parsed_variable = NORMALIZER_DEG;
            parsed_variables[NORMALIZER_DEG] = true;
          } else if (line == "normalize_to_den") {
            curr_parsed_variable = NORMALIZE_TO_DEN;
            parsed_variables[NORMALIZE_TO_DEN] = true;
          } else if (line == "normalizer_den_num") {
            curr_parsed_variable = NORMALIZER_DEN_NUM;
            parsed_variables[NORMALIZER_DEN_NUM] = true;
          } else if (line == "shifted_max_num_eqn") {
            curr_parsed_variable = SHIFTED_MAX_NUM_EQN;
            parsed_variables[SHIFTED_MAX_NUM_EQN] = true;
          } else if (line == "shift") {
            curr_parsed_variable = SHIFT;
            parsed_variables[SHIFT] = true;
          } else if (line == "shifted_degs_num") {
            curr_parsed_variable = SHIFTED_DEGS_NUM;
            parsed_variables[SHIFTED_DEGS_NUM] = true;
          } else if (line == "shifted_degs_den") {
            curr_parsed_variable = SHIFTED_DEGS_DEN;
            parsed_variables[SHIFTED_DEGS_DEN] = true;
          } else if (line == "zero_degs_num") {
            curr_parsed_variable = ZERO_DEGS_NUM;
            parsed_variables[ZERO_DEGS_NUM] = true;
          } else if (line == "zero_degs_den") {
            curr_parsed_variable = ZERO_DEGS_DEN;
            parsed_variables[ZERO_DEGS_DEN] = true;
          } else if (line == "g_ni") {
            curr_parsed_variable = G_NI;
            parsed_variables[G_NI] = true;
          } else if (line == "g_di") {
            curr_parsed_variable = G_DI;
            parsed_variables[G_DI] = true;
          } else if (line == "combined_ni") {
            curr_parsed_variable = COMBINED_NI;
            parsed_variables[COMBINED_NI] = true;
          } else if (line == "combined_di") {
            curr_parsed_variable = COMBINED_DI;
            parsed_variables[COMBINED_DI] = true;
          } else if (line == "interpolations") {
            curr_parsed_variable = INTERPOLATIONS;
            parsed_variables[INTERPOLATIONS] = true;
          } else {
            switch (curr_parsed_variable) {
              case COMBINED_PRIME: {
                combined_prime = mpz_class(line);
                break;
              }

              case IS_DONE: {
                std::unique_lock<std::mutex> lock_statics(mutex_status);
                done = std::stoi(line);
                break;
              }

              case MAX_DEG_NUM: {
                max_deg_num = std::stoi(line);
                break;
              }

              case MAX_DEG_DEN: {
                max_deg_den = std::stoi(line);
                break;
              }

              case NEED_PRIME_SHIFT: {
                if (!is_done())
                  need_prime_shift = std::stoi(line);

                break;
              }

              case NORMALIZER_DEG: {
                normalizer_deg = parse_vector(line);
                n = normalizer_deg.size();
                break;
              }

              case NORMALIZE_TO_DEN: {
                normalize_to_den = std::stoi(line);
                break;
              }

              case NORMALIZER_DEN_NUM: {
                normalizer_den_num = std::stoi(line);
                break;
              }

              case SHIFTED_MAX_NUM_EQN: {
                shifted_max_num_eqn = std::stoi(line);
                break;
              }

              case SHIFT: {
                if (!is_done()) {
                  std::vector<uint32_t> tmp_vec = parse_vector(line);
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);

                  shift = std::vector<FFInt> (n, 0);

                  for (uint32_t i = 0; i < n; ++i) {
                    if (tmp_vec[i] != 0)
                      shift[i] = FFInt(xorshift64star());
                  }
                }

                break;
              }

              case SHIFTED_DEGS_NUM: {
                std::vector<uint32_t> tmp_vec = parse_vector(line);

                for (const auto & el : tmp_vec) {
                  shifted_degs_num.emplace(el);
                }

                break;
              }

              case SHIFTED_DEGS_DEN: {
                std::vector<uint32_t> tmp_vec = parse_vector(line);

                for (const auto & el : tmp_vec) {
                  shifted_degs_den.emplace(el);
                }

                break;
              }

              case ZERO_DEGS_NUM: {
                std::vector<uint32_t> tmp_vec = parse_vector(line);

                for (const auto & el : tmp_vec) {
                  zero_degs_num.emplace(el);
                }

                break;
              }

              case ZERO_DEGS_DEN: {
                std::vector<uint32_t> tmp_vec = parse_vector(line);

                for (const auto & el : tmp_vec) {
                  zero_degs_den.emplace(el);
                }

                break;
              }

              case G_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                std::vector<mpz_class> tmp_rn = parse_rational_number(line);

                g_ni.emplace(std::make_pair(tmp_vec, RationalNumber(tmp_rn[0], tmp_rn[1])));

                break;
              }

              case G_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                std::vector<mpz_class> tmp_rn = parse_rational_number(line);

                g_di.emplace(std::make_pair(tmp_vec, RationalNumber(tmp_rn[0], tmp_rn[1])));
                break;

              }

              case COMBINED_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                combined_ni.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case COMBINED_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                combined_di.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case INTERPOLATIONS: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(-1);
                }

                interpolations = std::stoi(line);

                break;
              }
            }

          }
        }
      }

      for (const auto & el : parsed_variables) {
        if (!el) {
          ERROR_MSG("Incomplete input file! It cannot be used to resume a run.");
          std::exit(-1);
        }
      }

      file.close();

      for (const auto & el : combined_ni) add_non_solved_num(el.first);

      for (const auto & el : combined_di) add_non_solved_den(el.first);

      const_den = 0;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        is_singular_system = need_prime_shift;
      }

      std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
      new_prime = true;

      if (prime_number >= interpolations) {
        if (is_singular_system) {
          tmp_solved_coefs_den = 0;
          tmp_solved_coefs_num = 0;
          {
            std::unique_lock<std::mutex> lock_statics(mutex_statics);
            need_prime_shift = true;
          }

          for (const auto & el : g_ni) {
            if (normalize_to_den)
              add_non_solved_num(el.first);
            else if (el.first != std::vector<uint32_t> (n))
              add_non_solved_num(el.first);
          }

          for (const auto & el : g_di) {
            if (!normalize_to_den)
              add_non_solved_den(el.first);
            else if (el.first != std::vector<uint32_t> (n))
              add_non_solved_den(el.first);
          }

          if (normalize_to_den) {
            curr_deg_num = max_deg_num;

            if (max_deg_den > 0)
              curr_deg_den = max_deg_den;
          } else {
            curr_deg_den = max_deg_den;

            if (max_deg_num > 0)
              curr_deg_num = max_deg_num;
          }

          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = shifted_max_num_eqn;
        } else {
          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
        }
      }

      for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

      for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};
    } else {
      ERROR_MSG("The file '" + file_name + "' could not be found!");
      std::exit(-1);
    }
  }

  std::vector<uint32_t> RatReconst::parse_vector(std::string& line, int number_of_parameters) {
    size_t pos = 0;
    int i = 0;
    std::string delimiter = " ";
    std::vector<uint32_t> tmp {};

    if (number_of_parameters > 0)
      tmp.reserve(number_of_parameters);

    while ((pos = line.find(delimiter)) != std::string::npos) {
      tmp.emplace_back(std::stoi(line.substr(0, pos)));
      line.erase(0, pos + 1);
      i++;

      if (i == number_of_parameters) break;
    }

    return tmp;
  }

  std::vector<mpz_class> RatReconst::parse_rational_number(std::string& line) {
    size_t pos = line.find(" ");
    std::vector<mpz_class> tmp {};
    tmp.emplace_back(mpz_class(line.substr(0, pos)));
    line.erase(0, pos + 1);
    tmp.emplace_back(mpz_class(line));
    return tmp;
  }

  void RatReconst::parse_prime_number(std::string& file_name) {
    std::string reverse_file_name = file_name;
    std::reverse(reverse_file_name.begin(), reverse_file_name.end());
    reverse_file_name.erase(0, 4);
    size_t pos = reverse_file_name.find("_");
    prime_number = std::stoi(reverse_file_name.substr(0, pos)) + 1;
  }

  void RatReconst::set_singular_system_vars() {
    is_singular_system = true;
    tmp_solved_coefs_den = 0;
    tmp_solved_coefs_num = 0;

    for (const auto & el : g_ni) {
      if (normalize_to_den)
        add_non_solved_num(el.first);
      else if (el.first != std::vector<uint32_t> (n))
        add_non_solved_num(el.first);
    }

    for (const auto & el : g_di) {
      if (!normalize_to_den)
        add_non_solved_den(el.first);
      else if (el.first != std::vector<uint32_t> (n))
        add_non_solved_den(el.first);
    }

    if (normalize_to_den) {
      curr_deg_num = max_deg_num;

      if (max_deg_den > 0)
        curr_deg_den = max_deg_den;
    } else {
      curr_deg_den = max_deg_den;

      if (max_deg_num > 0)
        curr_deg_num = max_deg_num;
    }

    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = shifted_max_num_eqn;
    }

    for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

    for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};
  }

  void RatReconst::reset() {
    std::unique_lock<std::mutex> lock(mutex_statics);
    shift = std::vector<FFInt> ();
    need_prime_shift = false;
    set_singular_system = false;
    rand_zi = ff_pair_map();
    curr_shift = std::vector<uint32_t>();
    PolyReconst::reset();
  }

  void RatReconst::set_safe_interpolation() {
    interpolations = 100;
  }

  bool RatReconst::check_if_done(const FFInt& num, const FFInt& ti) {
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (!is_singular_system && set_singular_system) {
        lock_statics.unlock();
        set_singular_system_vars();
      }
    }

    if (rec_rat_coef()) {
      bool tmp_done = test_guess(num, ti);
      {
        std::unique_lock<std::mutex> lock(mutex_status);
        done = tmp_done;
      }

      if (done) {
        if (tag.size() > 0)
          save_state();

        std::unique_lock<std::mutex> lock(mutex_status);
        new_prime = false;
        curr_zi_order = std::vector<uint32_t>();
        use_chinese_remainder = false;
        return true;
      } else {
        for (const auto & ci : combined_ni) {
          g_ni.erase(ci.first);
        }

        for (const auto & ci : combined_di) {
          g_di.erase(ci.first);
        }
      }
    }

    if (!use_chinese_remainder) use_chinese_remainder = true;

    {
      std::unique_lock<std::mutex> lock(mutex_status);
      new_prime = false;
    }

    if (!is_singular_system) {
      ff_map tmp_num {};
      tmp_num.reserve(g_ni.size() + 1);
      ff_map tmp_den {};
      tmp_den.reserve(g_di.size() + 1);
      tmp_num.emplace(std::make_pair(std::vector<uint32_t> (n), 0));
      tmp_den.emplace(std::make_pair(std::vector<uint32_t> (n), 0));

      for (const auto & el : g_ni) {
        tmp_num[el.first] = FFInt(el.second.numerator) / FFInt(el.second.denominator);
      }

      for (const auto & el : g_di) {
        tmp_den[el.first] = FFInt(el.second.numerator) / FFInt(el.second.denominator);
      }

      solved_num = PolynomialFF(n, tmp_num);
      solved_den = PolynomialFF(n, tmp_den);
    }

    return false;
  }

  std::vector<FFInt> RatReconst::get_anchor_points() {
    std::vector<FFInt> res(n - 1);

    for (uint32_t i = 2; i <= n; ++i) {
      res[i - 2] = rand_zi[std::make_pair(i, 1)];
    }

    return res;
  }
}
