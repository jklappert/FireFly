//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "firefly/RatReconst.hpp"
#include "firefly/DenseSolver.hpp"
#include "firefly/gzstream.hpp"
#include "firefly/Logger.hpp"
#include "firefly/ParserUtils.hpp"
#include "firefly/ReconstHelper.hpp"
#include "firefly/utils.hpp"

#include <algorithm>
#include <sys/stat.h>

#ifdef FLINT
#include "/usr/include/flint/nmod_poly_factor.h"
#endif

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  std::unordered_set<uint32_t> RatReconst::singular_system_set {};
  std::vector<FFInt> RatReconst::global_anchor_points {};
  std::mutex RatReconst::mutex_statics;
  std::vector<uint32_t> RatReconst::curr_shift {};

  RatReconst::RatReconst(uint32_t n_) {
    n = n_;
    type = RAT;

    std::unique_lock<std::mutex> lock_status(mutex_status);

    combined_prime = mpz_class(std::to_string(FFInt::p));

    {
      std::lock_guard<std::mutex> lock_statics(mutex_statics);

      if (shift.empty()) {
        shift = std::vector<FFInt> (n);

        if (n > 1) {
          for (auto & el : shift) el = FFInt(get_rand_64());

          curr_shift = std::vector<uint32_t> (n, 1);
        }
      }
    }

    t_interpolator = ThieleInterpolator();

    if (n > 1) {
      curr_zi_order = std::vector<uint32_t> (n - 1, 1);
      lock_status.unlock();
      // add a zero to both polynomials to do arithmetics
      ff_map zero_deg {};
      zero_deg.emplace(std::make_pair(std::vector<uint32_t> (n), 0));
      solved_num = PolynomialFF(n, zero_deg);
      solved_den = PolynomialFF(n, zero_deg);

      // fill in the rand_vars for zi_order = 1
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (global_anchor_points.empty()) {
        lock_statics.unlock();
        generate_anchor_points();
      } else {
        private_anchor_points = global_anchor_points;
      }
    }
  }

  void RatReconst::scan_for_sparsest_shift() {
    scan = true;

    if (n == 1) {
      ERROR_MSG("You should never want to shift a univariate rational function.");
      std::exit(EXIT_FAILURE);
    }
  }

  void RatReconst::calc_factors(const std::string& var_, const std::pair<std::list<uint32_t>, std::list<uint32_t>>& degs) {
    is_calc_factors = true;
    var = var_;
    factor_degs = degs;

    if (!factor_degs.first.empty() || !factor_degs.second.empty())
      skip_thiele = true;

    for (const auto& el : factor_degs.first) {
      if (el == 0) {
        normalize_to_den = false;
        return;
      }
    }

    for (const auto& el : factor_degs.second) {
      if (el == 0) {
        normalize_to_den = true;
        return;
      }
    }
  }

  void RatReconst::set_zi_shift(const std::vector<uint32_t>& shifted_zis) {
    {
      std::lock_guard<std::mutex> lock(mutex_statics);
      shift = std::vector<FFInt> (n);
      curr_shift = shifted_zis;

      for (size_t i = 0; i != n; ++i) {
        if (shifted_zis[i] == 1)
          shift[i] = FFInt(get_rand_64());
      }
    }
  }

  std::pair<bool, bool> RatReconst::feed(const FFInt& new_ti, const FFInt& num,
                                         const std::vector<uint32_t>& fed_zi_ord,
                                         const uint32_t fed_prime) {
    std::lock_guard<std::mutex> lock(mutex_status);
    bool write_to_file = false;
    bool interpolate = false;

    if (!done && fed_prime == prime_number) {
      if (first_feed) {
        start_interpolation = true;
      }

      if (first_feed && !scan) {
        start = std::chrono::high_resolution_clock::now();

        if (num == 0) {
          new_prime = true;
          zero_counter ++;
          fed_zero = true;

          if (zero_counter == prime_number + 1)
            combined_prime = mpz_class(std::to_string(primes()[prime_number + 1]));

          if (tag.size() != 0) {
            if (prime_number != 0)
              save_zero_consecutive_prime();
            else
              save_zero_state();

            std::string probe_file_name_tmp = "ff_save/tmp/probe_" + tag + "_" + std::to_string(prime_number + 1) + "_tmp.gz";
	    std::rename(("ff_save/probes/" + tag + "_" + std::to_string(prime_number) + ".gz").c_str(), probe_file_name_tmp.c_str()); 
            ogzstream gzfile;
            std::string probe_file_name = "ff_save/probes/" + tag + "_" + std::to_string(prime_number + 1) + ".gz";
            gzfile.open(probe_file_name.c_str());
            gzfile.close();
	    std::rename(probe_file_name.c_str(), probe_file_name_tmp.c_str());
	    std::rename(probe_file_name_tmp.c_str(), probe_file_name.c_str());
          }

          ++prime_number;

          if (prime_number == 300) {
            ERROR_MSG("Your interpolation requests more than 300 primes.");
            std::exit(EXIT_FAILURE);
          } else if (zero_counter == 3 && prime_number == 3) {
            new_prime = false;
            done = true;
            g_ni[std::vector<uint32_t>(n)] = RationalNumber(0, 1);
            g_di[std::vector<uint32_t>(n)] = RationalNumber(1, 1);
          }
        } else {
          if (!check_interpolation)
            new_prime = false;

          first_feed = false;

          if (tag.size() != 0)
            saved_food.emplace(std::make_tuple(new_ti, num, fed_zi_ord));

          queue.emplace(std::make_tuple(new_ti, num, fed_zi_ord));
        }
      } else {
        if (!scan && tag.size() != 0)
          saved_food.emplace(std::make_tuple(new_ti, num, fed_zi_ord));

        queue.emplace(std::make_tuple(new_ti, num, fed_zi_ord));
      }

      if (start_interpolation) {
        start_interpolation = false;
        interpolate = true;
      }

      if (tag.size() != 0 && !scan && std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() > 1800.) {
        write_to_file = true;
      }
    }

    return std::make_pair(interpolate, write_to_file);
  }

  std::pair<bool, bool> RatReconst::feed(const std::vector<FFInt>& new_ti, const std::vector<FFInt>& num,
                                         const std::vector<std::vector<uint32_t>>& fed_zi_ord,
                                         const uint32_t fed_prime) {
    std::lock_guard<std::mutex> lock(mutex_status);
    bool write_to_file = false;
    bool interpolate = false;

    if (!done && fed_prime == prime_number) {
      if (first_feed) {
        start_interpolation = true;
      }

      if (first_feed && !scan) {
        start = std::chrono::high_resolution_clock::now();

        if (num[0] == 0) {
          new_prime = true;
          zero_counter ++;
          fed_zero = true;

          if (zero_counter == prime_number + 1)
            combined_prime = mpz_class(std::to_string(primes()[prime_number + 1]));

          if (tag.size() != 0) {
            if (prime_number != 0)
              save_zero_consecutive_prime();
            else
              save_zero_state();

            std::string probe_file_name_tmp = "ff_save/tmp/probe_" + tag + "_" + std::to_string(prime_number + 1) + "_tmp.gz";
	    std::rename(("ff_save/probes/" + tag + "_" + std::to_string(prime_number) + ".gz").c_str(), probe_file_name_tmp.c_str()); 
            ogzstream gzfile;
            std::string probe_file_name = "ff_save/probes/" + tag + "_" + std::to_string(prime_number + 1) + ".gz";
            gzfile.open(probe_file_name.c_str());
            gzfile.close();
	    std::rename(probe_file_name.c_str(), probe_file_name_tmp.c_str());
	    std::rename(probe_file_name_tmp.c_str(), probe_file_name.c_str());
          }

          ++prime_number;

          if (prime_number == 300) {
            ERROR_MSG("Your interpolation requests more than 300 primes.");
            std::exit(EXIT_FAILURE);
          } else if (zero_counter == 3 && prime_number == 3) {
            new_prime = false;
            done = true;
            g_ni[std::vector<uint32_t>(n)] = RationalNumber(0, 1);
            g_di[std::vector<uint32_t>(n)] = RationalNumber(1, 1);
          }
        } else {
          if (!check_interpolation)
            new_prime = false;

          first_feed = false;

          for (size_t i = 0; i != new_ti.size(); ++i) {
            if (tag.size() != 0)
              saved_food.emplace(std::make_tuple(new_ti[i], num[i], fed_zi_ord[i]));

            queue.emplace(std::make_tuple(new_ti[i], num[i], fed_zi_ord[i]));
          }
        }
      } else {
        for (size_t i = 0; i != new_ti.size(); ++i) {
          if (!scan && tag.size() != 0)
            saved_food.emplace(std::make_tuple(new_ti[i], num[i], fed_zi_ord[i]));

          queue.emplace(std::make_tuple(new_ti[i], num[i], fed_zi_ord[i]));
        }
      }

      if (start_interpolation) {
        start_interpolation = false;
        interpolate = true;
      }

      if (tag.size() != 0 && !scan && std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() > 1800.) {
        write_to_file = true;
      }
    }

    return std::make_pair(interpolate, write_to_file);
  }

  // TODO Change return type since the first bool should not be needed anymore
  std::tuple<bool, bool, uint32_t> RatReconst::interpolate() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (fed_zero) {
      fed_zero = false;
    } else if (is_interpolating || queue.empty()) {
      return std::make_tuple(false, done, prime_number);
    } else {
      is_interpolating = true;

      while (!queue.empty()) {
        auto food = queue.front();
        queue.pop();

        lock.unlock();

        interpolate(std::get<0>(food), std::get<1>(food), std::get<2>(food));

        while (saved_ti.find(curr_zi_order) != saved_ti.end()) {
          lock.lock();

          std::pair<FFInt, FFInt> key_val = saved_ti[curr_zi_order].front();
          saved_ti[curr_zi_order].pop();

          lock.unlock();

          interpolate(key_val.first, key_val.second, curr_zi_order);
        }

        lock.lock();
      }

      is_interpolating = false;
      start_interpolation = true;
    }

    return std::make_tuple(true, done, prime_number);
  }

  void RatReconst::interpolate(const FFInt& new_ti,
                               const FFInt& num,
                               const std::vector<uint32_t>& fed_zi_ord) {
    if (!done) {
      // Compare if the food is the expected food; if not, store it for later use
      if (fed_zi_ord == curr_zi_order) {
        // first check if we are done. If not start the interpolation again using
        // the chinese remainder theorem in combining the previous results
        if (new_prime)
          if (check_if_done(num, new_ti))
            return;

        if (max_deg_num == -1) {// Use Thiele
          if (!skip_thiele)
            check = t_interpolator.add_point(num, new_ti);
          else
            check = build_factor_uni_gauss(new_ti, num);
        } else {
          std::queue<std::pair<FFInt, FFInt>> t_food;

          // Prepare food for Gauss system
          if (n > 1) {
            if (saved_ti.find(curr_zi_order) != saved_ti.end()) {
              std::lock_guard<std::mutex> lock(mutex_status);

              t_food = saved_ti[curr_zi_order];
              saved_ti.erase(curr_zi_order);
            }
          }

          t_food.emplace(std::make_pair(new_ti, num));

          // Iterate through all feeds and build the uni/multivariate Gauss
          // system
          while (!t_food.empty()) {
            auto food = t_food.front();
            t_food.pop();
            FFInt tmp_ti = food.first;
            FFInt tmp_num = food.second;

            // Get yi's for the current feed
            std::vector<FFInt> yis;

            if (n > 1)
              yis = get_rand_zi_vec(curr_zi_order);

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
            if (!skip_thiele)
              canonical = t_interpolator.get_result();
            else
              canonical = solve_factor_uni_gauss();

            PolynomialFF numerator = PolynomialFF(1, canonical.first);
            PolynomialFF denominator = PolynomialFF(1, canonical.second);

            if (scan) {
              {
                std::lock_guard<std::mutex> lock(mutex_status);
                done = true;
              }

              uint32_t tmp_max_deg_num = numerator.max_deg()[0];
              uint32_t tmp_max_deg_den = denominator.max_deg()[0];
              factor_degs.first.clear();
              factor_degs.second.clear();

              for (const auto& el : numerator.coefs) {
                factor_degs.first.emplace_back(el.first[0]);
              }

              for (const auto& el : denominator.coefs) {
                factor_degs.second.emplace_back(el.first[0]);
              }

              {
                std::lock_guard<std::mutex> lock_status(mutex_status);
                shift_works = false;
              }

              if (curr_shift == std::vector<uint32_t> (n, 1)) {
                all_shift_max_degs = {tmp_max_deg_num, tmp_max_deg_den};
              } else if (tmp_max_deg_num == all_shift_max_degs[0] && tmp_max_deg_den == all_shift_max_degs[1]) {
                if (denominator.min_deg()[0] == 0)
                  normalize_to_den = true;
                else
                  normalize_to_den = false;

                shift_works = true;
              }

              t_interpolator = ThieleInterpolator();
              // TODO return list of degrees for system of equations
              // similar to get result ff, clear when shift is accepted no external input needed
              return;
            }

            factor_degs.first.clear();
            factor_degs.second.clear();

            max_deg_num = numerator.max_deg()[0];
            max_deg_den = denominator.max_deg()[0];

            if (n > 1 && shift != std::vector<FFInt>(n)) {
              PolynomialFF zero_poly(n, {{std::vector<uint32_t> (n), 0}});

              // Initialize subtraction terms with zero
              for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_num); ++i) {
                sub_num[i] = zero_poly;
              }

              for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_den); ++i) {
                sub_den[i] = zero_poly;
              }
            }

            if (n == 1) {
              normalizer_deg = denominator.min_deg();
              normalizer_den_num = true;
            }

            FFInt equalizer;

            if (normalize_to_den)
              equalizer = FFInt(1) / denominator.coefs[denominator.min_deg()];
            else
              equalizer = FFInt(1) / numerator.coefs[numerator.min_deg()];

            canonical.first = (numerator * equalizer).coefs;
            canonical.second = (denominator * equalizer).coefs;

            if (n > 1) {
              // Remove constants and vanishing degrees from system of equations
              for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_num); ++i) {
                if (canonical.first.find( {i}) == canonical.first.end())
                  tmp_solved_coefs_num ++;
                else if (i == 0 && canonical.first.find( {0}) != canonical.first.end()) {
                  tmp_solved_coefs_num ++;
                  std::vector<uint32_t> zero_deg(1);
                  solved_degs_num[0] = PolynomialFF(n, {{std::vector<uint32_t> (n), canonical.first[zero_deg]}});
                  canonical.first.erase(zero_deg);
                  dense_solve_degs_num.emplace(0);
                }
                else
                  dense_solve_degs_num.emplace(i); // First write all in dense and remove if required
              }

              for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_den); ++i) {
                if (canonical.second.find( {i}) == canonical.second.end())
                  tmp_solved_coefs_den ++;
                else if (i == 0 && canonical.second.find( {0}) != canonical.second.end()) {
                  tmp_solved_coefs_den ++;
                  std::vector<uint32_t> zero_deg(1);
                  solved_degs_den[0] = PolynomialFF(n, {{std::vector<uint32_t> (n), canonical.second[zero_deg]}});
                  canonical.second.erase(zero_deg);
                  dense_solve_degs_den.emplace(0);
                }
                else
                  dense_solve_degs_den.emplace(i); // First write all in dense and remove if required
              }

              // set number of equations needed for univariate rational function
              // reconstruction needed for multivariate polynomial feed
              {
                std::lock_guard<std::mutex> lock(mutex_status);
                num_eqn = max_deg_den + max_deg_num + 2 - tmp_solved_coefs_num - tmp_solved_coefs_den;
              }
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
            // remove zero degrees
            if (first_run) {
              ff_map num_coef = canonical.first;
              ff_map den_coef = canonical.second;
              std::vector<uint32_t> degs;

              for (const auto & el : num_coef) {
                if (el.second == 0) {
                  canonical.first.erase(el.first);
                  tmp_solved_coefs_num ++;
                } else
                  degs.emplace_back(el.first[0]);
              }

              std::sort(degs.begin(), degs.end(),
              [](const uint32_t & l, const uint32_t & r) {
                return l > r;
              });

              for (const auto & deg : degs) {
                coef_n.emplace_back(std::make_pair(deg, PolyReconst(n - 1, deg, true)));

                if (is_set_individual_degree_bounds) {
                  coef_n.back().second.set_individual_degree_bounds(get_individual_degree_bounds(false));
                }
              }

              degs.clear();

              for (const auto & el : den_coef) {
                if (el.second == 0) {
                  canonical.second.erase(el.first);
                  tmp_solved_coefs_den ++;
                } else
                  degs.emplace_back(el.first[0]);
              }

              std::sort(degs.begin(), degs.end(),
              [](const uint32_t & l, const uint32_t & r) {
                return l > r;
              });

              for (const auto & deg : degs) {
                coef_d.emplace_back(std::make_pair(deg, PolyReconst(n - 1, deg, true)));

                if (is_set_individual_degree_bounds) {
                  coef_d.back().second.set_individual_degree_bounds(get_individual_degree_bounds(true));
                }
              }

              std::lock_guard<std::mutex> lock(mutex_status);
              zi = 2;
            }

            // Save the current results to the map to access them later
            // Numerator
            for (const auto & el : canonical.first) {
              uint32_t deg = el.first[0];

              if (curr_zi_order[zi - 2] < deg + 2)
                saved_num_num[curr_zi_order][ {deg, zi}] = el.second;
            }

            // Denominator
            for (const auto & el : canonical.second) {
              uint32_t deg = el.first[0];

              if (curr_zi_order[zi - 2] < deg + 2)
                saved_num_den[curr_zi_order][ {deg, zi}] = el.second;
            }

            if (first_run) first_run = false;

            // Interpolate the numerator. Always try to interpolate all degrees.
            // Prefer higher zi and zi_order to remove as many degrees as possible
            // from the system of equations. When a higher degree is interpolated,
            // try to interpolate the missing degrees sparsely.
            for (auto it = coef_n.begin(); it != coef_n.end(); ++it) {
              if (it->second.get_zi_order() == curr_zi_order || restart_sparse_interpolation_num) {
                bool poly_done = false;

                if (restart_sparse_interpolation_num) {
                  it -> second = PolyReconst(n - 1, it->first, true);

                  if (is_set_individual_degree_bounds) {
                    it->second.set_individual_degree_bounds(get_individual_degree_bounds(false));
                  }
                  restart_sparse_interpolation_num = false;
                }

                if (dense_solve_degs_num.find(it->first) != dense_solve_degs_num.end())
                  poly_done = feed_poly(it->first, it->second, saved_num_num, true);
                else
                  poly_done = feed_poly(it->first, it->second, saved_num_num, true, true);

                if (poly_done) {
                  it = coef_n.erase(it);
                  --it;
                }
              }
            }

            for (auto it = coef_d.begin(); it != coef_d.end(); ++it) {
              if (it->second.get_zi_order() == curr_zi_order || restart_sparse_interpolation_den) {
                bool poly_done = false;

                if (restart_sparse_interpolation_den) {
                  it -> second = PolyReconst(n - 1, it->first, true);

                  if (is_set_individual_degree_bounds) {
                    it->second.set_individual_degree_bounds(get_individual_degree_bounds(true));
                  }
                  restart_sparse_interpolation_den = false;
                }

                if (dense_solve_degs_den.find(it->first) != dense_solve_degs_den.end())
                  poly_done = feed_poly(it->first, it->second, saved_num_den, false);
                else
                  poly_done = feed_poly(it->first, it->second, saved_num_den, false, true);

                if (poly_done) {
                  it = coef_d.erase(it);
                  --it;
                }
              }
            }

            // First reset values to give a starting point for comparison
            {
              std::lock_guard<std::mutex> lock(mutex_status);
              zi = 2;
              curr_zi_order = std::vector<uint32_t> (n - 1, 1);
              vandermonde_size = 1;
            }

            // After all new interpolation points, iterate over all interpolation
            // objects and check which should be feeded next
            // Numerator
            for (auto & p_rec : coef_n) {
              uint32_t tmp_zi = p_rec.second.get_zi() + 1;
              std::vector<uint32_t> tmp_zi_ord = p_rec.second.get_zi_order();

              if (zi < tmp_zi) {
                std::lock_guard<std::mutex> lock(mutex_status);
                zi = tmp_zi;
                curr_zi_order = tmp_zi_ord;
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              } else if (zi == tmp_zi && a_grt_b(tmp_zi_ord, curr_zi_order)) {
                std::lock_guard<std::mutex> lock(mutex_status);
                curr_zi_order = tmp_zi_ord;
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              } else if (zi == tmp_zi && a_eq_b(tmp_zi_ord, curr_zi_order) && vandermonde_size < p_rec.second.get_vandermonde_num_eqn()) {
                std::unique_lock<std::mutex> lock(mutex_status);
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              }
            }

            // Denominator
            for (auto & p_rec : coef_d) {
              uint32_t tmp_zi = p_rec.second.get_zi() + 1;
              std::vector<uint32_t> tmp_zi_ord = p_rec.second.get_zi_order();

              if (zi < tmp_zi) {
                std::lock_guard<std::mutex> lock(mutex_status);
                zi = tmp_zi;
                curr_zi_order = tmp_zi_ord;
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              } else if (zi == tmp_zi && a_grt_b(tmp_zi_ord, curr_zi_order)) {
                std::lock_guard<std::mutex> lock(mutex_status);
                curr_zi_order = tmp_zi_ord;
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              } else if (zi == tmp_zi && a_eq_b(tmp_zi_ord, curr_zi_order) && vandermonde_size < p_rec.second.get_vandermonde_num_eqn()) {
                std::unique_lock<std::mutex> lock(mutex_status);
                vandermonde_size = p_rec.second.get_vandermonde_num_eqn();
              }
            }

            // combine results
            if (coef_n.size() == 0 && coef_d.size() == 0) {
              saved_num_num.clear();
              saved_num_den.clear();
              vandermonde_size = 1;
              FFInt const_den = 0;

              // Calculate shift polynomials and combine with previous ones if there is any shift
              // To do so, we have to remove the shift from all solved degrees
              if (shift != std::vector<FFInt> (n)) {
                for (int i = max_deg_num; i > -1; i--) {
                  if (dense_solve_degs_num.find(static_cast<uint32_t>(i)) != dense_solve_degs_num.end()) {
                    solved_degs_num[i] -= sub_num[i];

                    if (i != 0) {
                      for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_num[i], i)) {
                        sub_num[tmp_shift.first] += tmp_shift.second;
                      }
                    }
                  }
                }

                for (int i = max_deg_den; i > -1; i--) {
                  if (dense_solve_degs_den.find(static_cast<uint32_t>(i)) != dense_solve_degs_den.end()) {
                    solved_degs_den[i] -= sub_den[i];

                    if (i != 0) {
                      for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_den[i], i)) {
                        sub_den[tmp_shift.first] += tmp_shift.second;
                      }
                    }
                  }
                }

                if (normalize_to_den)
                  const_den = sub_den[0].coefs[std::vector<uint32_t> (n)];
                else
                  const_den = sub_num[0].coefs[std::vector<uint32_t> (n)];
              }

              first_run = true;

              // Remove normalization due to the shift
              PolynomialFF numerator;
              PolynomialFF denominator;

              for (auto & el : solved_degs_num) {
                PolynomialFF res = el.second;

                // TODO: define empty/zero polynomials uniquely
                if (!(res.coefs.size() == 0) && !(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  numerator += res;
                else
                  zero_degs_num.emplace(el.first);
              }

              for (auto & el : solved_degs_den) {
                PolynomialFF res = el.second;

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
              if (prime_number == 0 + zero_counter) {
                if (const_den != 1) {
                  terminator = FFInt(1) - const_den;

                  if (normalize_to_den)
                    normalizer_den_num = true;
                  else
                    normalizer_den_num = false;

                  normalizer_deg = std::vector<uint32_t> (n, 0);
                } else if (!normalize_to_den && numerator.coefs.find(std::vector<uint32_t> (n, 0)) != numerator.coefs.end()) {
                  terminator = numerator.coefs[std::vector<uint32_t> (n, 0)];
                  normalizer_den_num = false;
                  normalizer_deg = std::vector<uint32_t> (n, 0);
                } else if (normalize_to_den && denominator.coefs.find(std::vector<uint32_t> (n, 0)) != denominator.coefs.end()) {
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

                // -1 for the normalization constant
                shifted_max_num_eqn = solved_degs_num.size() + solved_degs_den.size() - 1;
              } else if (normalizer_den_num)
                terminator = denominator.coefs[normalizer_deg];
              else
                terminator = numerator.coefs[normalizer_deg];

              dense_solve_degs_den.clear();
              dense_solve_degs_num.clear();
              solved_degs_num.clear();
              solved_degs_den.clear();

              if (interpolations > prime_number + 1)
                is_singular_system = false;

              // normalize
              FFInt equalizer = FFInt(1) / terminator;

              if (equalizer == 0)
                div_by_zero = true;

              numerator *= equalizer;
              denominator *= equalizer;

              combine_primes(numerator.coefs, denominator.coefs);

              std::lock_guard<std::mutex> lock(mutex_status);
              std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
              zi = 2;
            }

            return;
          } else {
            // build multivariate Vandermonde systems and evaluate them if possible
            if (!is_singular_system) {
              for (const auto & sol : canonical.first) {
                uint32_t key = sol.first[0];

                if (coef_mat_num[key].size() < non_solved_degs_num[key].size())
                  coef_mat_num[key].emplace_back(sol.second);

                // Solve multivariate Vandermonde system for corresponding degree,
                // remove entry from non_solved_degs and add it to solve_degs
                if (coef_mat_num[key].size() == non_solved_degs_num[key].size()) {
                  if (solved_degs_num[key].zero())
                    solved_degs_num[key] = solve_vandermonde_system(non_solved_degs_num[key], coef_mat_num[key], get_anchor_points());
                  else
                    solved_degs_num[key] += solve_vandermonde_system(non_solved_degs_num[key], coef_mat_num[key], get_anchor_points());

                  non_solved_degs_num.erase(key);
                  coef_mat_num.erase(key);
                }
              }

              for (const auto & sol : canonical.second) {
                uint32_t key = sol.first[0];

                if (coef_mat_den[key].size() < non_solved_degs_den[key].size())
                  coef_mat_den[key].emplace_back(sol.second);

                if (coef_mat_den[key].size() == non_solved_degs_den[key].size()) {
                  if (solved_degs_den[key].zero())
                    solved_degs_den[key] = solve_vandermonde_system(non_solved_degs_den[key], coef_mat_den[key], get_anchor_points());
                  else
                    solved_degs_den[key] += solve_vandermonde_system(non_solved_degs_den[key], coef_mat_den[key], get_anchor_points());

                  non_solved_degs_den.erase(key);
                  coef_mat_den.erase(key);
                }
              }
            } else {
              for (const auto & sol : canonical.first) {
                uint32_t key = sol.first[0];

                if (coef_mat_num[key].size() < non_solved_degs_num[key].size())
                  coef_mat_num[key].emplace_back(sol.second);
              }

              // Need to iterate two times trhough this map to ensure that any equation has all terms
              for (const auto & sol : canonical.first) {
                uint32_t key = sol.first[0];

                // Solve multivariate Vandermonde system for corresponding degree,
                // remove entry from non_solved_degs and add it to solve_degs
                if ((dense_solve_degs_num.find(key) != dense_solve_degs_num.end() || key == max_num_coef_num.first) && coef_mat_num[key].size() == non_solved_degs_num[key].size()) {
                  // When we reach this point, we can solve all remaining polynomials sparsely

                  if (key == max_num_coef_num.first) {
                    for (uint32_t i = max_deg_num; i > 0; i--) {
                      if (coef_mat_num.find(i) != coef_mat_num.end()) {
                        // Check if this is a zero polynomial
                        if (zero_degs_num.find(i) != zero_degs_num.end()) {
                          solved_degs_num[i] = PolynomialFF(n, {{std::vector<uint32_t> (n), 0}});
                          non_solved_degs_num.erase(i);
                          coef_mat_num.erase(i);
                          tmp_solved_coefs_num ++;
                          continue;
                        }

                        // Remove shift
                        if (i < static_cast<uint32_t>(max_deg_num)) {
                          for (uint32_t j = 0; j < coef_mat_num[i].size(); ++j) {
                            std::vector<uint32_t> tmp_zi_ord(n - 1 , j + 1);
                            coef_mat_num[i][j] -= sub_num[i].calc_n_m_1(get_rand_zi_vec(tmp_zi_ord));
                          }
                        }

                        solved_degs_num[i] = solve_vandermonde_system(non_solved_degs_num[i], coef_mat_num[i], get_anchor_points());
                        non_solved_degs_num.erase(i);
                        coef_mat_num.erase(i);
                        tmp_solved_coefs_num ++;

                        PolynomialFF tmp_poly = solved_degs_num[i];

                        // Calculate shift and combine with previous
                        for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_num[i], i)) {
                          sub_num[tmp_shift.first] += tmp_shift.second;
                        }
                      } else if (i > 0 && solved_degs_num.find(i) != solved_degs_num.end() && zero_degs_num.find(i) == zero_degs_num.end()) { // when already solved we just need to calculate the shift
                        PolynomialFF tmp_poly = solved_degs_num[i];
                        tmp_poly -= sub_num[i];

                        // Calculate shift and combine with previous
                        for (const auto & tmp_shift : calculate_shift_polynomials(tmp_poly, i)) {
                          sub_num[tmp_shift.first] += tmp_shift.second;
                        }
                      } else if (i != static_cast<uint32_t>(max_deg_num) && !sub_num[i].zero() && solved_degs_num.find(i) == solved_degs_num.end()) {
                        solved_degs_num[i] = sub_num[i];
                        non_solved_degs_num.erase(i);
                        tmp_solved_coefs_num ++;
                      }
                    }

                    if (coef_mat_num.find(0) != coef_mat_num.end()) {
                      // Check if this is a zero polynomial
                      if (zero_degs_num.find(0) != zero_degs_num.end()) {
                        solved_degs_num[0] = PolynomialFF(n, {{std::vector<uint32_t> (n), 0}});
                        non_solved_degs_num.erase(0);
                        coef_mat_num.erase(0);
                        tmp_solved_coefs_num ++;
                      } else {

                        // Remove shift
                        if (0 < max_deg_num) {
                          std::vector<FFInt> zero_vec(n - 1 , 0);
                          coef_mat_num[0][0] -= sub_num[0].calc_n_m_1(zero_vec);
                        }

                        solved_degs_num[0] = solve_vandermonde_system(non_solved_degs_num[0], coef_mat_num[0], get_anchor_points());
                        non_solved_degs_num.erase(0);
                        coef_mat_num.erase(0);
                        tmp_solved_coefs_num ++;
                      }
                    }

                    break;
                  }

                  solved_degs_num[key] = solve_vandermonde_system(non_solved_degs_num[key], coef_mat_num[key], get_anchor_points());

                  non_solved_degs_num.erase(key);
                  coef_mat_num.erase(key);
                  tmp_solved_coefs_num ++;
                }
              }

              for (const auto & sol : canonical.second) {
                uint32_t key = sol.first[0];

                if (coef_mat_den[key].size() < non_solved_degs_den[key].size())
                  coef_mat_den[key].emplace_back(sol.second);
              }

              for (const auto & sol : canonical.second) {
                uint32_t key = sol.first[0];

                if ((dense_solve_degs_den.find(key) != dense_solve_degs_den.end() || key == max_num_coef_den.first) && coef_mat_den[key].size() == non_solved_degs_den[key].size()) {

                  // When we reach this point, we can solve all remaining polynomials sparsely
                  if (key == max_num_coef_den.first) {
                    for (uint32_t i = max_deg_den; i > 0; i--) {
                      if (coef_mat_den.find(i) != coef_mat_den.end()) {
                        // Check if this is a zero polynomial
                        if (zero_degs_den.find(i) != zero_degs_den.end()) {
                          solved_degs_den[i] = PolynomialFF(n, {{std::vector<uint32_t> (n), 0}});
                          non_solved_degs_den.erase(i);
                          coef_mat_den.erase(i);
                          tmp_solved_coefs_den ++;
                          continue;
                        }

                        // Remove shift
                        if (i < static_cast<uint32_t>(max_deg_den)) {
                          for (uint32_t j = 0; j < coef_mat_den[i].size(); ++j) {
                            std::vector<uint32_t> tmp_zi_ord(n - 1 , j + 1);
                            coef_mat_den[i][j] -= sub_den[i].calc_n_m_1(get_rand_zi_vec(tmp_zi_ord));
                          }
                        }

                        solved_degs_den[i] = solve_vandermonde_system(non_solved_degs_den[i], coef_mat_den[i], get_anchor_points());
                        non_solved_degs_den.erase(i);
                        coef_mat_den.erase(i);
                        tmp_solved_coefs_den ++;

                        // Calculate shift and combine with previous
                        for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_den[i], i)) {
                          sub_den[tmp_shift.first] += tmp_shift.second;
                        }
                      } else if (i > 0 && solved_degs_den.find(i) != solved_degs_den.end() && zero_degs_den.find(i) == zero_degs_den.end()) { // when already solved we just need to calculate the shift
                        PolynomialFF tmp_poly = solved_degs_den[i];
                        tmp_poly -= sub_den[i];

                        // Calculate shift and combine with previous
                        for (const auto & tmp_shift : calculate_shift_polynomials(tmp_poly, i)) {
                          sub_den[tmp_shift.first] += tmp_shift.second;
                        }
                      } else if (i != static_cast<uint32_t>(max_deg_den) && !sub_den[i].zero() && solved_degs_den.find(i) == solved_degs_den.end()) {
                        solved_degs_den[i] = sub_den[i];
                        non_solved_degs_den.erase(i);
                        tmp_solved_coefs_den ++;
                      }

                      if (coef_mat_den.find(0) != coef_mat_den.end()) {
                        // Check if this is a zero polynomial
                        if (zero_degs_den.find(0) != zero_degs_den.end()) {
                          solved_degs_den[0] = PolynomialFF(n, {{std::vector<uint32_t> (n), 0}});
                          non_solved_degs_den.erase(0);
                          coef_mat_den.erase(0);
                          tmp_solved_coefs_den ++;
                        } else {

                          // Remove shift
                          if (0 < max_deg_den) {
                            std::vector<FFInt> zero_vec(n - 1 , 0);
                            coef_mat_den[0][0] -= sub_den[0].calc_n_m_1(zero_vec);
                          }

                          solved_degs_den[0] = solve_vandermonde_system(non_solved_degs_den[0], coef_mat_den[0], get_anchor_points());
                          non_solved_degs_den.erase(0);
                          coef_mat_den.erase(0);
                          tmp_solved_coefs_den ++;
                        }
                      }
                    }

                    break;
                  }

                  solved_degs_den[key] = solve_vandermonde_system(non_solved_degs_den[key], coef_mat_den[key], get_anchor_points());
                  non_solved_degs_den.erase(key);
                  coef_mat_den.erase(key);
                  tmp_solved_coefs_den ++;
                }
              }
            }

            // promote to next prime and combine results
            if (coef_mat_num.empty() && coef_mat_den.empty()) {
              if (is_singular_system) {
                non_solved_degs_den.clear();
                non_solved_degs_num.clear();

                for (const auto & el : solved_degs_num) {
                  // Only add polynomials which are not zero when the shift is subtracted
                  if (zero_degs_num.find(el.first) == zero_degs_num.end()) {
                    if (dense_solve_degs_num.find(el.first) != dense_solve_degs_num.end()) {
                      PolynomialFF tmp_poly = el.second - sub_num[el.first];
                      solved_num += tmp_poly;

                      for (const auto & coef : tmp_poly.coefs) {
                        if (coef.second != 0)
                          non_solved_degs_num[el.first].emplace_back(coef.first);
                      }
                    } else {
                      solved_num += el.second;

                      for (const auto & coef : el.second.coefs) {
                        if (coef.second != 0)

                          non_solved_degs_num[el.first].emplace_back(coef.first);
                      }
                    }
                  }
                }

                for (const auto & el : solved_degs_den) {
                  // Only add polynomials which are not zero when the shift is subtracted
                  if (zero_degs_den.find(el.first) == zero_degs_den.end()) {
                    if (dense_solve_degs_den.find(el.first) != dense_solve_degs_den.end()) {
                      PolynomialFF tmp_poly = el.second - sub_den[el.first];
                      solved_den += tmp_poly;

                      for (const auto & coef : tmp_poly.coefs) {
                        if (coef.second != 0)
                          non_solved_degs_den[el.first].emplace_back(coef.first);
                      }
                    } else {
                      solved_den += el.second;

                      for (const auto & coef : el.second.coefs) {
                        if (coef.second != 0)
                          non_solved_degs_den[el.first].emplace_back(coef.first);
                      }
                    }
                  }
                }

                // normalize
                FFInt terminator = 0;

                FFInt const_den = 0;

                if (normalize_to_den)
                  const_den = sub_den[0].coefs[std::vector<uint32_t> (n)];
                else
                  const_den = sub_num[0].coefs[std::vector<uint32_t> (n)];

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
              } else {
                for (const auto & el : solved_degs_num) {
                  solved_num += el.second;
                }

                for (const auto & el : solved_degs_den) {
                  solved_den += el.second;
                }
              }

              // remove the constant if it is zero
              //solved_num.remove_zero_coefs();
              //solved_den.remove_zero_coefs();
              if (solved_num.coefs.find(std::vector<uint32_t> (n, 0)) != solved_num.coefs.end() && solved_num.coefs[std::vector<uint32_t> (n, 0)] == 0) {
                solved_num.coefs.erase(std::vector<uint32_t> (n, 0));
              }

              if (solved_den.coefs.find(std::vector<uint32_t> (n, 0)) != solved_den.coefs.end() && solved_den.coefs[std::vector<uint32_t> (n, 0)] == 0) {
                solved_den.coefs.erase(std::vector<uint32_t> (n, 0));
              }

              combine_primes(solved_num.coefs, solved_den.coefs);
              {
                std::lock_guard<std::mutex> lock(mutex_status);
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
              std::lock_guard<std::mutex> lock(mutex_status);
              saved_ti.erase(curr_zi_order);
              std::transform(curr_zi_order.begin(), curr_zi_order.end(),
              curr_zi_order.begin(), [](uint32_t x) {return x + 1;});

              if (!is_singular_system)
                num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
              else
                num_eqn = shifted_max_num_eqn - tmp_solved_coefs_num - tmp_solved_coefs_den;
            }
          }
        }
      } else {
        if (saved_ti.find(fed_zi_ord) == saved_ti.end()) {
          std::lock_guard<std::mutex> lock(mutex_status);

          std::queue<std::pair<FFInt, FFInt>> tmp_;
          tmp_.emplace(std::make_pair(new_ti, num));
          saved_ti[fed_zi_ord] = tmp_;
        } else {
          std::lock_guard<std::mutex> lock(mutex_status);

          saved_ti[fed_zi_ord].emplace(std::make_pair(new_ti, num));
        }
      }
    }
  }

  bool RatReconst::feed_poly(uint32_t curr_deg,
                             PolyReconst& rec,
                             ff_map_map& saved_num,
                             bool is_num,
                             bool sparse
                            ) {
    uint32_t tmp_zi = rec.get_zi() + 1;
    std::vector<uint32_t> tmp_zi_ord = rec.get_zi_order();

    while (!rec.is_new_prime()) {
      try {
        std::vector<uint32_t> key = {curr_deg, tmp_zi};
        FFInt food = saved_num.at(tmp_zi_ord).at(key);
        // delete unused saved data
        //saved_num[tmp_zi_ord].erase(key);
        // set random values for the yis
        std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);

        // feed to PolyReconst
        if (sparse && shift != std::vector<FFInt> (n)) {
          FFInt sub = 0;

          if (is_num)
            sub = sub_num[curr_deg].calc_n_m_1(yis);
          else
            sub = sub_den[curr_deg].calc_n_m_1(yis);

          rec.feed(food - sub);
        } else
          rec.feed(food);

        tmp_zi = rec.get_zi() + 1;
        tmp_zi_ord = rec.get_zi_order();
      } catch (std::out_of_range& e) {
        return false;
      }
    }

    if (is_num) {
      solved_degs_num[curr_deg] = rec.get_result_ff();
      tmp_solved_coefs_num ++;

      if (shift != std::vector<FFInt> (n)) {
        if (curr_deg == coef_n.front().first) {
          for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_num[curr_deg], curr_deg)) {
            sub_num[tmp_shift.first] += tmp_shift.second;
          }

          if (dense_solve_degs_num.find(curr_deg) != dense_solve_degs_num.end())
            dense_solve_degs_num.erase(curr_deg);

          // Reset the polynomial interplolation and do it sparsely
          if (coef_n.size() > 1) {
            uint32_t tmp_deg = (++coef_n.begin()) -> first;
            dense_solve_degs_num.erase(tmp_deg);
            restart_sparse_interpolation_num = true;

            for (uint32_t tmp_deg_2 = curr_deg - 1; tmp_deg_2 > tmp_deg; tmp_deg_2--) {
              if (solved_degs_num.find(tmp_deg_2) != solved_degs_num.end()) {
                solved_degs_num[tmp_deg_2] -= sub_num[tmp_deg_2];

                if (dense_solve_degs_num.find(tmp_deg_2) != dense_solve_degs_num.end())
                  dense_solve_degs_num.erase(tmp_deg_2);

                for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_num[tmp_deg_2], tmp_deg_2)) {
                  sub_num[tmp_shift.first] += tmp_shift.second;
                }
              }
            }
          }
        }
      }
    } else {
      solved_degs_den[curr_deg] = rec.get_result_ff();
      tmp_solved_coefs_den ++;

      if (shift != std::vector<FFInt> (n)) {
        if (curr_deg == coef_d.front().first) {
          for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_den[curr_deg], curr_deg)) {
            sub_den[tmp_shift.first] += tmp_shift.second;
          }

          if (dense_solve_degs_den.find(curr_deg) != dense_solve_degs_den.end())
            dense_solve_degs_den.erase(curr_deg);

          // Reset the polynomial interplolation and do it sparsely
          if (coef_d.size() > 1) {
            uint32_t tmp_deg = (++coef_d.begin()) -> first;
            dense_solve_degs_den.erase(tmp_deg);
            restart_sparse_interpolation_den = true;

            for (uint32_t tmp_deg_2 = curr_deg - 1; tmp_deg_2 > tmp_deg; tmp_deg_2--) {
              if (solved_degs_den.find(tmp_deg_2) != solved_degs_den.end()) {

                solved_degs_den[tmp_deg_2] -= sub_den[tmp_deg_2];

                if (dense_solve_degs_den.find(tmp_deg_2) != dense_solve_degs_den.end())
                  dense_solve_degs_den.erase(tmp_deg_2);

                for (const auto & tmp_shift : calculate_shift_polynomials(solved_degs_den[tmp_deg_2], tmp_deg_2)) {
                  sub_den[tmp_shift.first] += tmp_shift.second;
                }
              }
            }
          }
        }
      }
    }

    {
      std::lock_guard<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 2 - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }

    return true;
  }

  void RatReconst::combine_primes(ff_map& numerator, ff_map& denominator) {
    std::vector<uint32_t> tmp_deg_num {};
    std::vector<uint32_t> tmp_deg_den {};

    // store temporary results
    tmp_res_num = numerator;
    tmp_res_den = denominator;

    {
      std::lock_guard<std::mutex> lock(mutex_status);

      saved_ti = ff_queue_map();
    }

    max_num_coef_num = std::make_pair(0, 0);
    max_num_coef_den = std::make_pair(0, 0);

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
    saved_num_num = ff_map_map();
    saved_num_den = ff_map_map();
    bool is_safe_mode = interpolations != 1;

    if (!div_by_zero) {
      if (!use_chinese_remainder) {
        combined_prime = mpz_class(std::to_string(FFInt::p));
        combined_ni = convert_to_mpz(numerator);
        combined_di = convert_to_mpz(denominator);

        mpz_map combined_ni_back = combined_ni;

        for (auto & c_ni : combined_ni_back) {

          if (!normalizer_den_num && c_ni.first == normalizer_deg)
            remove_ni(c_ni.first, RationalNumber(1, 1));
          else {
            add_non_solved_num(c_ni.first);

            if (is_safe_mode)
              combined_primes_ni[c_ni.first] = combined_prime;
          }
        }

        mpz_map combined_di_back = combined_di;

        for (auto & c_di : combined_di_back) {

          if (normalizer_den_num && c_di.first == normalizer_deg)
            remove_di(c_di.first, RationalNumber(1, 1));
          else {
            add_non_solved_den(c_di.first);

            if (is_safe_mode)
              combined_primes_di[c_di.first] = combined_prime;
          }
        }

        if (is_singular_system) {
          check_for_solved_degs(tmp_deg_num, true);

          if (is_singular_system)
            check_for_solved_degs(tmp_deg_den, false);
        }
      } else {
        mpz_map combined_ni_back = combined_ni;
        mpz_map combined_di_back = combined_di;

        if (is_safe_mode) {
          for (const auto & el : numerator) {
            if (combined_ni.find(el.first) == combined_ni.end() && g_ni.find(el.first) == g_ni.end())
              combined_ni.emplace(std::make_pair(el.first, 0));
          }

          for (const auto & el : denominator) {
            if (combined_di.find(el.first) == combined_di.end() && g_di.find(el.first) == g_di.end())
              combined_di.emplace(std::make_pair(el.first, 0));
          }
        }

        mpz_class combined_prime_back = combined_prime;
        mpz_map combined_prime_ni_back = combined_primes_ni;
        mpz_map combined_prime_di_back = combined_primes_di;
        std::pair<mpz_class, mpz_class> p1;
        std::pair<mpz_class, mpz_class> p2;
        std::pair<mpz_class, mpz_class> p3;

        //numerator
        for (auto it = combined_ni.begin(); it != combined_ni.end(); ++it) {
          if (!is_safe_mode) {
            p1 = std::make_pair(it->second, combined_prime);
            p2 = std::make_pair(mpz_class(std::to_string(numerator[it->first].n)), mpz_class(std::to_string(FFInt::p)));
            p3 = run_chinese_remainder(p1, p2);
            combined_ni[it->first] = p3.first;
          } else {
            if (numerator.find(it->first) != numerator.end() && g_ni.find(it->first) == g_ni.end()) {
              if (combined_ni[it->first] == 0) {
                combined_ni[it->first] = mpz_class(std::to_string(numerator[it->first].n));
                combined_primes_ni[it->first] = mpz_class(std::to_string(FFInt::p));
              } else {
                p1 = std::make_pair(it->second, combined_primes_ni[it->first]);
                p2 = std::make_pair(mpz_class(std::to_string(numerator[it->first].n)), mpz_class(std::to_string(FFInt::p)));
                p3 = run_chinese_remainder(p1, p2);
                combined_ni[it->first] = p3.first;
                combined_primes_ni[it->first] = p3.second;
              }
            }
          }
        }

        // denominator
        for (auto it = combined_di.begin(); it != combined_di.end(); ++it) {
          if (!is_safe_mode) {
            p1 = std::make_pair(it->second, combined_prime);
            p2 = std::make_pair(mpz_class(std::to_string(denominator[it->first].n)), mpz_class(std::to_string(FFInt::p)));
            p3 = run_chinese_remainder(p1, p2);
            combined_di[it->first] = p3.first;
          } else {
            if (denominator.find(it->first) != denominator.end() && g_di.find(it->first) == g_di.end()) {
              if (combined_di[it->first] == 0) {
                combined_di[it->first] = mpz_class(std::to_string(denominator[it->first].n));
                combined_primes_di[it->first] = mpz_class(std::to_string(FFInt::p));
              } else {
                p1 = std::make_pair(it->second, combined_primes_di[it->first]);
                p2 = std::make_pair(mpz_class(std::to_string(denominator[it->first].n)), mpz_class(std::to_string(FFInt::p)));
                p3 = run_chinese_remainder(p1, p2);
                combined_di[it->first] = p3.first;
                combined_primes_di[it->first] = p3.second;
              }
            }
          }
        }

        combined_prime = p3.second;

        // Remove already known coefficients from solve algorithm to save numerical runs
        for (auto & c_ni : combined_ni_back) {
          if (g_ni.find(c_ni.first) == g_ni.end() && (!is_safe_mode || (is_safe_mode && numerator.find(c_ni.first) != numerator.end()))) {
            RationalNumber last_rn_wang;
            RationalNumber curr_rn_wang;
            bool last_wang;
            bool curr_wang;

            std::pair<bool, RationalNumber> res;

            if (!is_safe_mode)
              res = get_rational_coef(c_ni.second, combined_prime_back);
            else
              res = get_rational_coef(c_ni.second, combined_prime_ni_back[c_ni.first]);

            last_wang = res.first;

            if (last_wang)
              last_rn_wang = res.second;

            if (!is_safe_mode)
              res = get_rational_coef(combined_ni[c_ni.first], combined_prime);
            else
              res = get_rational_coef(combined_ni[c_ni.first], combined_primes_ni[c_ni.first]);

            curr_wang = res.first;

            if (curr_wang)
              curr_rn_wang = res.second;

            RationalNumber last_rn_monagan;
            RationalNumber curr_rn_monagan;
            bool last_monagan;
            bool curr_monagan;

            if (!is_safe_mode)
              res = get_rational_coef_mqrr(c_ni.second, combined_prime_back);
            else
              res = get_rational_coef_mqrr(c_ni.second, combined_prime_ni_back[c_ni.first]);

            last_monagan = res.first;

            if (last_monagan)
              last_rn_monagan = res.second;

            if (!is_safe_mode)
              res = get_rational_coef_mqrr(combined_ni[c_ni.first], combined_prime);
            else
              res = get_rational_coef_mqrr(combined_ni[c_ni.first], combined_primes_ni[c_ni.first]);

            curr_monagan = res.first;

            if (curr_monagan)
              curr_rn_monagan = res.second;

            if (last_wang && curr_wang && last_rn_wang == curr_rn_wang) {
              remove_ni(c_ni.first, curr_rn_wang);

              if (is_safe_mode)
                combined_primes_ni.erase(c_ni.first);

              continue;
            }

            if (last_monagan && curr_monagan && last_rn_monagan == curr_rn_monagan) {
              remove_ni(c_ni.first, curr_rn_monagan);

              if (is_safe_mode)
                combined_primes_ni.erase(c_ni.first);

              continue;
            }

            if (!is_safe_mode && c_ni.second == combined_ni[c_ni.first]) {
              RationalNumber rn = RationalNumber(c_ni.second, 1);
              remove_ni(c_ni.first, rn);
            } else
              add_non_solved_num(c_ni.first);
          }
        }

        for (auto & c_di : combined_di_back) {
          if (g_di.find(c_di.first) == g_di.end() && (!is_safe_mode || (is_safe_mode && denominator.find(c_di.first) != denominator.end()))) {
            RationalNumber last_rn_wang;
            RationalNumber curr_rn_wang;
            bool last_wang;
            bool curr_wang;
            std::pair<bool, RationalNumber> res;

            if (!is_safe_mode)
              res = get_rational_coef(c_di.second, combined_prime_back);
            else
              res = get_rational_coef(c_di.second, combined_prime_di_back[c_di.first]);

            last_wang = res.first;

            if (last_wang)
              last_rn_wang = res.second;

            if (!is_safe_mode)
              res = get_rational_coef(combined_di[c_di.first], combined_prime);
            else
              res = get_rational_coef(combined_di[c_di.first], combined_primes_di[c_di.first]);

            curr_wang = res.first;

            if (curr_wang)
              curr_rn_wang = res.second;


            RationalNumber last_rn_monagan;
            RationalNumber curr_rn_monagan;
            bool last_monagan;
            bool curr_monagan;

            if (!is_safe_mode)
              res = get_rational_coef_mqrr(c_di.second, combined_prime_back);
            else
              res = get_rational_coef_mqrr(c_di.second, combined_prime_di_back[c_di.first]);

            last_monagan = res.first;

            if (last_monagan)
              last_rn_monagan = res.second;

            if (!is_safe_mode)
              res = get_rational_coef_mqrr(combined_di[c_di.first], combined_prime);
            else
              res = get_rational_coef_mqrr(combined_di[c_di.first], combined_primes_di[c_di.first]);

            curr_monagan = res.first;

            if (curr_monagan)
              curr_rn_monagan = res.second;

            if (last_wang && curr_wang && last_rn_wang == curr_rn_wang) {
              remove_di(c_di.first, curr_rn_wang);

              if (is_safe_mode)
                combined_primes_di.erase(c_di.first);

              continue;
            }

            if (last_monagan && curr_monagan && last_rn_monagan == curr_rn_monagan) {
              remove_di(c_di.first, curr_rn_monagan);

              if (is_safe_mode)
                combined_primes_di.erase(c_di.first);

              continue;
            }

            if (!is_safe_mode && c_di.second == combined_di[c_di.first]) {
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
        // Check if the state should be written out after this prime
        if (tag.size() != 0)
          save_state();

        if (is_singular_system) {
          set_singular_system_vars();
          {
            std::lock_guard<std::mutex> lock_statics(mutex_statics);

            if (singular_system_set.find(prime_number + 1) == singular_system_set.end())
              singular_system_set.emplace(prime_number + 1);
          }

          std::lock_guard<std::mutex> lock(mutex_status);
          num_eqn = shifted_max_num_eqn;
        } else {
          if (n > 1) {
            for (const auto & el : non_solved_degs_num) {
              coef_mat_num[el.first] = std::vector<FFInt> {};
            }

            for (const auto & el : non_solved_degs_den) {
              coef_mat_den[el.first] = std::vector<FFInt> {};
            }
          }

          std::lock_guard<std::mutex> lock(mutex_status);
          num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
        }
      }
    } else {
      div_by_zero = false;
      interpolations ++;
    }

    if (prime_number + 1 < interpolations) {
      max_deg_num = -1;
      first_run = true;

      // Check if the state should be written out after this prime
      if (tag.size() != 0)
        save_state();
    }

    {
      std::unique_lock<std::mutex> lock(mutex_status);

      ++prime_number;

      // Remove old probes and create new file
      if (tag.size() != 0) {
        while (is_writing_probes) {
          // TODO: better with condition variable: as pointer?
          lock.unlock();
          lock.lock();
        }

        std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> tmp_;
        saved_food.swap(tmp_);

	std::string probe_file_name_tmp = "ff_save/tmp/probe_" + tag + "_" + std::to_string(prime_number) + "_tmp.gz";
	std::rename(("ff_save/probes/" + tag + "_" + std::to_string(prime_number - 1) + ".gz").c_str(), probe_file_name_tmp.c_str()); 
	ogzstream gzfile;
	std::string probe_file_name = "ff_save/probes/" + tag + "_" + std::to_string(prime_number) + ".gz";
	gzfile.open(probe_file_name.c_str());
	gzfile.close();
	std::rename(probe_file_name.c_str(), probe_file_name_tmp.c_str());
	std::rename(probe_file_name_tmp.c_str(), probe_file_name.c_str());
      }

      queue = std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>>();
      saved_ti.clear();
      zi = 1;
      new_prime = true;
      first_feed = true;
      check_interpolation = true;

      if (prime_number == 300) {
        ERROR_MSG("Your interpolation requests more than 300 primes.");
        std::exit(EXIT_FAILURE);
      }
    }

    tmp_solved_coefs_den = 0;
    tmp_solved_coefs_num = 0;

#ifdef FLINT
    if (is_calc_factors) {
        nmod_poly_t numerator, denominator;
        nmod_poly_factor_t fac_numerator, fac_denominator;
        nmod_poly_init(numerator, FFInt::p);
        nmod_poly_init(denominator, FFInt::p);
        nmod_poly_factor_init(fac_numerator);
        nmod_poly_factor_init(fac_denominator);

        // Rewrite numerator in FLINTs notation and determine maximal degree
        for (const auto & coef : get_result_ff().numerator.coefs) {
          if (coef.first[0] > static_cast<uint32_t>(max_deg_num))
            max_deg_num = coef.first[0];

          nmod_poly_set_coeff_ui(numerator, coef.first[0], coef.second.n);
        }

        // Rewrite denominator in FLINTs notation and determine maximal degree
        for (const auto & coef : get_result_ff().denominator.coefs) {
          if (coef.first[0] > static_cast<uint32_t>(max_deg_den))
            max_deg_den = coef.first[0];

          nmod_poly_set_coeff_ui(denominator, coef.first[0], coef.second.n);
        }

        // Factorize polynomials
        //TODO replace FLINT with NTL? According to NTLs benchmarks their factorization algorithms are more efficient
        nmod_poly_factor_kaltofen_shoup(fac_numerator, numerator);
        nmod_poly_factor_kaltofen_shoup(fac_denominator, denominator);
        //nmod_poly_factor_cantor_zassenhaus(fac_numerator, numerator);
        //nmod_poly_factor_cantor_zassenhaus(fac_denominator, denominator);
        // Berlekamp is fastest but does not find all factors for some reason
        //nmod_poly_factor_berlekamp(fac_numerator, numerator);
        //nmod_poly_factor_berlekamp(fac_denominator, denominator);

        // Rewrite and store in result objects
        if (fac_numerator[0].num + fac_denominator[0].num != 0) {
          factors.first.reserve(fac_numerator[0].num);
          factors_flint.first.reserve(fac_numerator[0].num);
          for (int fac_num_num = 0; fac_num_num != fac_numerator[0].num; ++fac_num_num) {
            char *tmp_fac_str_pretty_ptr = nmod_poly_get_str_pretty(&fac_numerator[0].p[fac_num_num], var.c_str());
            std::string tmp_fac(tmp_fac_str_pretty_ptr);
            std::string tmp_factor = "(" + tmp_fac + ")^" + std::to_string(fac_numerator[0].exp[fac_num_num]);
            factors.first.emplace_back(tmp_factor);
            flint_free(tmp_fac_str_pretty_ptr);
            char *tmp_fac_str_ptr = nmod_poly_get_str(&fac_numerator[0].p[fac_num_num]);
            std::string tmp_fac_str(tmp_fac_str_ptr);
            factors_flint.first.emplace_back(std::make_pair(tmp_fac_str, fac_numerator[0].exp[fac_num_num]));
            flint_free(tmp_fac_str_ptr);
          }

          factors.second.reserve(fac_denominator[0].num);
          factors_flint.second.reserve(fac_denominator[0].num);
          for (int fac_den_num = 0; fac_den_num != fac_denominator[0].num; ++fac_den_num) {
            char *tmp_fac_str_pretty_ptr = nmod_poly_get_str_pretty(&fac_denominator[0].p[fac_den_num], var.c_str());
            std::string tmp_fac(tmp_fac_str_pretty_ptr);
            std::string tmp_factor = "(" + tmp_fac + ")^" + std::to_string(fac_denominator[0].exp[fac_den_num]);
            factors.second.emplace_back(tmp_factor);
            flint_free(tmp_fac_str_pretty_ptr);
            char *tmp_fac_str_ptr = nmod_poly_get_str(&fac_denominator[0].p[fac_den_num]);
            std::string tmp_fac_str(tmp_fac_str_ptr);
            factors_flint.second.emplace_back(std::make_pair(tmp_fac_str, fac_denominator[0].exp[fac_den_num]));
            flint_free(tmp_fac_str_ptr);
          }
        }

        // Free memory
        nmod_poly_clear(numerator);
        nmod_poly_clear(denominator);
        nmod_poly_factor_clear(fac_numerator);
        nmod_poly_factor_clear(fac_denominator);
    }
#endif
  }

  RationalFunction RatReconst::get_result() {
    std::lock_guard<std::mutex> lock(mutex_status);

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
      std::exit(EXIT_FAILURE);
    }
  }

  RationalFunctionFF RatReconst::get_result_ff() {
    std::lock_guard<std::mutex> lock(mutex_status);

    ff_map tmp_num = convert_to_ffint(g_ni);
    ff_map tmp_den = convert_to_ffint(g_di);

    tmp_res_num.insert(tmp_num.begin(), tmp_num.end());
    tmp_res_den.insert(tmp_den.begin(), tmp_den.end());

    return RationalFunctionFF(PolynomialFF(n, tmp_res_num), PolynomialFF(n, tmp_res_den));
  }

#ifdef FLINT
  std::pair<std::vector<std::string>, std::vector<std::string>> RatReconst::get_factors_ff() {
    return factors;
  }

  // TODO Workaround to be not called from main thread
  std::pair<std::unordered_map<uint32_t, uint64_t>, std::unordered_map<uint32_t, uint64_t>> RatReconst::get_canonical_factors(const std::pair<std::unordered_set<uint32_t>, std::unordered_set<uint32_t>>& factors_pos) {
    nmod_poly_t c_numerator, c_denominator;
    nmod_poly_init(c_numerator, FFInt::p);
    nmod_poly_init(c_denominator, FFInt::p);
    nmod_poly_set_coeff_ui(c_numerator, 0, 1);
    nmod_poly_set_coeff_ui(c_denominator, 0, 1);

    // Numerator
    for (size_t i = 0; i != factors_flint.first.size(); ++i) {
      if (factors_pos.first.find(i) != factors_pos.first.end()) {
        for (size_t j = 0; j != factors_flint.first[i].second; ++j) {
          nmod_poly_t tmp_1, tmp_2;
          nmod_poly_init(tmp_1, FFInt::p);
          nmod_poly_init(tmp_2, FFInt::p);
          nmod_poly_swap(tmp_1, c_numerator);
          nmod_poly_set_str(tmp_2, factors_flint.first[i].first.c_str());
          nmod_poly_mul_classical(c_numerator, tmp_1, tmp_2);
          nmod_poly_clear(tmp_1);
          nmod_poly_clear(tmp_2);
        }
      }
    }

    // Denominator
    for (size_t i = 0; i != factors_flint.second.size(); ++i) {
      if (factors_pos.second.find(i) != factors_pos.second.end()) {
        for (size_t j = 0; j != factors_flint.second[i].second; ++j) {
          nmod_poly_t tmp_1, tmp_2;
          nmod_poly_init(tmp_1, FFInt::p);
          nmod_poly_init(tmp_2, FFInt::p);
          nmod_poly_swap(tmp_1, c_denominator);
          nmod_poly_set_str(tmp_2, factors_flint.second[i].first.c_str());
          nmod_poly_mul_classical(c_denominator, tmp_1, tmp_2);
          nmod_poly_clear(tmp_1);
          nmod_poly_clear(tmp_2);
        }
      }
    }

    // Set result
    std::unordered_map<uint32_t, uint64_t> c_numerator_map {};
    std::unordered_map<uint32_t, uint64_t> c_denominator_map {};

    // Numerator
    if (factors_pos.first.size() != 0) {
      for (size_t i = 0; i != static_cast<size_t>(c_numerator->length); ++i) {
        if (c_numerator->coeffs[i])
          c_numerator_map.emplace(std::make_pair(i, c_numerator->coeffs[i]));
      }
    }

    // Denominator
    if (factors_pos.second.size() != 0) {
      for (size_t i = 0; i != static_cast<size_t>(c_denominator->length); ++i) {
        if (c_denominator->coeffs[i])
          c_denominator_map.emplace(std::make_pair(i, c_denominator->coeffs[i]));
      }
    }

    // Free memory
    nmod_poly_clear(c_numerator);
    nmod_poly_clear(c_denominator);
    return std::make_pair(c_numerator_map, c_denominator_map);
  }
#endif

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
      std::vector<uint32_t> power = {el.first};
      numerator.emplace(std::make_pair(std::move(power), results[counter]));
      counter ++;
    }

    for (const auto & el : coef_d) {
      std::vector<uint32_t> power = {el.first};
      denominator.emplace(std::make_pair(std::move(power), results[counter]));
      counter ++;
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
      std::lock_guard<std::mutex> lock_statics(mutex_statics);
      yis[0] = ti + shift[0];

      for (uint32_t i = 1; i != n; ++i) {
        yis[i] = ti * get_rand_zi(i + 1, curr_zi_order[i - 1]) + shift[i];
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
    std::lock_guard<std::mutex> lock_statics(mutex_statics);
    shift = std::vector<FFInt> (n, 0);
  }

  bool RatReconst::build_factor_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num) {
    std::vector<FFInt> eq;
    eq.reserve(factor_degs.first.size() + factor_degs.second.size());

    FFInt res(0);
    // Build system of equations;
    for (const auto & el : factor_degs.first) {
      if (!scan && !normalize_to_den && el == 0)
        res = -1;
      else if (scan && el == all_shift_max_degs[0])
        res = -tmp_ti.pow(el);
      else
        eq.emplace_back(tmp_ti.pow(el));
    }

    for (const auto & el : factor_degs.second) {
      if (!scan && normalize_to_den && el == 0)
        res = tmp_num;
      else
        eq.emplace_back(-tmp_num * tmp_ti.pow(el));
    }

    eq.emplace_back(res);
    coef_mat.emplace_back(std::move(eq));

    size_t size = factor_degs.first.size() + factor_degs.second.size() - 1;

    if (coef_mat.size() == size)
      return true;
    else
      return false;
  }

  std::pair<ff_map, ff_map> RatReconst::solve_factor_uni_gauss() {
    std::vector<FFInt> results = solve_gauss_system(coef_mat, coef_mat.size());
    coef_mat.clear();

    ff_map numerator;
    ff_map denominator;
    uint32_t counter = 0;

    for (const auto & el : factor_degs.first) {
      std::vector<uint32_t> power = {el};
      if ((!scan && !normalize_to_den && el == 0) || (scan && el == all_shift_max_degs[0]))
        numerator.emplace(std::make_pair(std::move(power), 1));
      else {
        if (results[counter] != 0)
          numerator.emplace(std::make_pair(std::move(power), results[counter]));

        counter ++;
      }
    }

    for (const auto & el : factor_degs.second) {
      std::vector<uint32_t> power = {el};
      if (!scan && normalize_to_den && el == 0)
        denominator.emplace(std::make_pair(std::move(power), 1));
      else {
        if (results[counter] != 0)
          denominator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }
    }

    return std::make_pair(numerator, denominator);
  }

  void RatReconst::build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);

    // Build system of equations;
    for (const auto & el : coef_n) {
      eq.emplace_back(tmp_ti.pow(el.first));
    }

    for (const auto & el : coef_d) {
      eq.emplace_back(-tmp_num * tmp_ti.pow(el.first));
    }

    FFInt res(0);

    // Build result vector including subtracted coefficients which have already
    // been solved
    if (coef_mat.size() == 0) {
      for (auto & el : solved_degs_num) {
        uint32_t tmp_deg = el.first;

        if (dense_solve_degs_num.find(tmp_deg) != dense_solve_degs_num.end())
          num_sub_num[tmp_deg] = el.second.calc_n_m_1(yis);
        else
          num_sub_num[tmp_deg] = el.second.calc_n_m_1(yis) + sub_num[tmp_deg].calc_n_m_1(yis);
      }

      for (auto & el : solved_degs_den) {
        uint32_t tmp_deg = el.first;

        if (dense_solve_degs_den.find(tmp_deg) != dense_solve_degs_den.end())
          num_sub_den[tmp_deg] = el.second.calc_n_m_1(yis);
        else
          num_sub_den[tmp_deg] = el.second.calc_n_m_1(yis) + sub_den[tmp_deg].calc_n_m_1(yis);
      }
    }

    for (const auto & el : num_sub_num) {
      res -= el.second * tmp_ti.pow(el.first);
    }

    for (const auto & el : num_sub_den) {
      res += el.second * tmp_ti.pow(el.first) * tmp_num;
    }

    eq.emplace_back(res);
    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);
    FFInt res(0);

    // Build system of equations; in combined_.. are the non-solved coefficients
    for (const auto & pow_vec : non_solved_degs_num) {
      eq.emplace_back(tmp_ti.pow(pow_vec.first));
    }

    for (const auto & pow_vec : non_solved_degs_den) {
      eq.emplace_back(-tmp_num * tmp_ti.pow(pow_vec.first));
    }

    if (!is_singular_system) {
      // Build result vector including subtracted coefficients which have already
      // been solved
      if (coef_mat.size() == 0) {
        for (auto & el : solved_degs_num) {
          num_sub_num[el.first] = el.second.calc_n_m_1(yis);
        }

        for (auto & el : solved_degs_den) {
          num_sub_den[el.first] = el.second.calc_n_m_1(yis);
        }
      }
    } else {
      if (normalize_to_den)
        res = tmp_num;
      else
        res = -1;

      // Build result vector including subtracted coefficients which have already
      // been solved
      if (coef_mat.size() == 0) {
        for (auto & el : solved_degs_num) {
          uint32_t tmp_deg = el.first;

          if (dense_solve_degs_num.find(tmp_deg) != dense_solve_degs_num.end())
            num_sub_num[tmp_deg] = el.second.calc_n_m_1(yis);
          else
            num_sub_num[tmp_deg] = el.second.calc_n_m_1(yis) + sub_num[tmp_deg].calc_n_m_1(yis);
        }

        for (auto & el : solved_degs_den) {
          uint32_t tmp_deg = el.first;

          if (dense_solve_degs_den.find(tmp_deg) != dense_solve_degs_den.end())
            num_sub_den[tmp_deg] = el.second.calc_n_m_1(yis);
          else
            num_sub_den[tmp_deg] = el.second.calc_n_m_1(yis) + sub_den[tmp_deg].calc_n_m_1(yis);
        }
      }
    }

    for (const auto & el : num_sub_num) {
      res -= el.second * tmp_ti.pow(el.first);
    }

    for (const auto & el : num_sub_den) {
      res += el.second * tmp_ti.pow(el.first) * tmp_num;
    }

    eq.emplace_back(res);
    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::generate_anchor_points() {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);

    global_anchor_points = std::vector<FFInt> (n - 1, 0);

    for (uint32_t tmp_zi = 0; tmp_zi != n - 1; ++tmp_zi) {
      global_anchor_points[tmp_zi] = FFInt(get_rand_64());
    }

    private_anchor_points = global_anchor_points;

    PolyReconst rec(n - 1, 0, true);
    rec.set_anchor_points(private_anchor_points, true);

    if (!shift.empty() && n > 1) {
      for (auto & el : shift) el = FFInt(el.n);
    }
  }

  void RatReconst::set_anchor_points(const std::vector<FFInt>& anchor_points) {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);

    global_anchor_points = anchor_points;
    private_anchor_points = anchor_points;

    PolyReconst rec(n - 1, 0, true);
    rec.set_anchor_points(anchor_points, true);
  }

  void RatReconst::set_shift(const std::vector<FFInt>& shift_) {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);

    if (n > 1)
      shift = shift_;
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
    if (max_deg_num == 0 || max_deg_den == 0) {
      is_singular_system = false;
      return;
    }

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

  FFInt RatReconst::get_rand_zi(uint32_t zi, uint32_t order) const {
    return private_anchor_points[zi - 2].pow(order);
  }

  std::vector<FFInt> RatReconst::get_rand_zi_vec(const std::vector<uint32_t>& order) {
    std::vector<FFInt> res {};

    for (size_t i = 0; i != n - 1; ++i) {
      res.emplace_back(private_anchor_points[i].pow(order[i]));
    }

    return res;
  }

  FFInt RatReconst::get_zi_shift(uint32_t zi) const {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);
    return shift[zi - 1];
  }

  bool RatReconst::is_shift_working() {
    std::lock_guard<std::mutex> lock_status(mutex_status);
    std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
    done = false;
    zi = 1;
    return shift_works;
  }

  void RatReconst::accept_shift() {
    if (tag.size() != 0) {
      std::string file_name = "ff_save/states/" + tag + "_" + std::to_string(prime_number) + ".gz";
      //std::ofstream file;
      ogzstream file;
      file.open(file_name.c_str());
      file << "tag_name\n" << tag_name << "\n";
      file << "normalize_to_den\n" << std::to_string(normalize_to_den) << "\n";
      file << "individual_degrees_num\n";

      if (!individual_degree_bounds.empty()) {
        for (const auto & el : individual_degree_bounds) {
          file << el.first << " ";
        }

        file << "\n";
      }

      file << "individual_degrees_den\n";

      if (!individual_degree_bounds.empty()) {
        for (const auto & el : individual_degree_bounds) {
          file << el.second << " ";
        }

        file << "\n";
      }

      file.close();
    }

    /*factor_degs.first.clear();
    factor_degs.second.clear();*/
    skip_thiele = true;

    scan = false;
  }

  bool RatReconst::get_is_interpolating() const {
    std::lock_guard<std::mutex> lock(mutex_status);
    return is_interpolating;
  }

  std::vector<FFInt> RatReconst::get_zi_shift_vec() const {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);
    return shift;
  }

  bool RatReconst::need_shift(uint32_t prime_counter) {
    std::lock_guard<std::mutex> lock_statics(mutex_statics);

    for (uint32_t i = 0; i < n; ++i) {
      shift[i] = FFInt(shift[i].n);
    }

    if (singular_system_set.find(prime_counter) == singular_system_set.end())
      return false;
    else
      return true;
  }

  polff_map RatReconst::calculate_shift_polynomials(const PolynomialFF& poly, uint32_t deg) {
    polff_map res;
    std::vector<uint32_t> zero_deg(n);
    PolynomialFF zero_poly(n, {{zero_deg, 0}});

    for (uint32_t i = 0; i != deg; ++i) {
      res.emplace(std::make_pair(i, zero_poly));
    }

    std::vector<FFInt> tmp_shift;
    {
      std::lock_guard<std::mutex> lock_statics(mutex_statics);
      tmp_shift = shift;
    }

    PolynomialFF sub_pol = poly.add_shift(tmp_shift);

    for (auto & el : sub_pol.coefs) {
      uint32_t tmp_deg = 0;

      for (const auto & deg_ : el.first) tmp_deg += deg_;

      if (tmp_deg < deg) {
        if (el.second > 0)
          res[tmp_deg] += PolynomialFF(n, {{el.first, el.second}});
      }
    }

    return res;
  }

  void RatReconst::set_tag(const std::string& tag_) {
    if (tag.size() == 0) {
      tag = tag_;
    } else
      WARNING_MSG("This object has already a valid tag!");
  }

  void RatReconst::set_tag_name(const std::string& tag_name_) {
    if (tag_name.size() == 0) {
      tag_name = tag_name_;
      std::string file_name = "ff_save/states/" + tag + "_" + std::to_string(prime_number) + ".gz";
      ogzstream file;
      file.open(file_name.c_str());
      file << "tag_name\n" << tag_name << "\n";
      file << "normalize_to_den\n1\n";
      file << "individual_degrees_num\n";

      if (!individual_degree_bounds.empty()) {
        for (const auto & el : individual_degree_bounds) {
          file << el.first << " ";
        }

        file << "\n";
      }

      file << "individual_degrees_den\n";

      if (!individual_degree_bounds.empty()) {
        for (const auto & el : individual_degree_bounds) {
          file << el.second << " ";
        }

        file << "\n";
      }

      file.close();
      ogzstream file_2;
      file_name = "ff_save/probes/" + tag + "_" + std::to_string(prime_number) + ".gz";
      file_2.open(file_name.c_str());
      file_2.close();
    } else
      WARNING_MSG("This object has already a valid tag!");
  }

  void RatReconst::save_zero_state() {
    ogzstream file;
    std::string file_name = "ff_save/states/" + tag + "_" + std::to_string(prime_number) + ".gz";
    file.open(file_name.c_str());
    file << "ZERO\n";
    interpolations > 1 ? file << "1\n" : file << "0\n";
    file << tag_name << "\n";
    file.close();
  }

  void RatReconst::save_zero_consecutive_prime() {
    std::string file_name_old = "ff_save/states/" + tag + "_" + std::to_string(prime_number - 1) + ".gz";
    std::string file_name_new = "ff_save/states/" + tag + "_" + std::to_string(prime_number) + ".gz";

    if (std::rename(file_name_old.c_str(), file_name_new.c_str()) != 0)
      WARNING_MSG("The previously saved file '" + file_name_old + "' could not be renamed.");
  }

  void RatReconst::save_state() {
    ogzstream file;
    std::string file_name = "ff_save/tmp/" + tag + "_" + std::to_string(prime_number) + ".gz";
    file.open(file_name.c_str());
    file << "combined_prime\n" << combined_prime.get_str() << "\n";
    file << "tag_name\n" << tag_name << "\n";
    file << "is_done\n" << is_done() << "\n";
    file << "max_deg_num\n" << max_deg_num << "\n";
    file << "max_deg_den\n" << max_deg_den << "\n";
    file << "individual_degrees_num\n";

    if (!individual_degree_bounds.empty()) {
      for (const auto & el : individual_degree_bounds) {
        file << el.first << " ";
      }

      file << "\n";
    }

    file << "individual_degrees_den\n";

    if (!individual_degree_bounds.empty()) {
      for (const auto & el : individual_degree_bounds) {
        file << el.second << " ";
      }

      file << "\n";
    }

    if (interpolations != 1)
      file << "need_prime_shift\n" << "1" << "\n";
    else
      file << "need_prime_shift\n" << is_singular_system << "\n";

    file << "normalizer_deg\n";

    for (uint32_t i = 0; i != n; ++i) {
      file << normalizer_deg[i] << " ";
    }

    file << "\n";

    file << "normalize_to_den\n" << normalize_to_den << "\n";
    file << "normalizer_den_num\n" << normalizer_den_num << "\n";
    file << "shifted_max_num_eqn\n" << shifted_max_num_eqn << "\n";
    file << "shift\n";

    {
      std::lock_guard<std::mutex> lock(mutex_statics);

      for (uint32_t i = 0; i != n; ++i) {
        file << (shift[i] != 0 ? "1 " : "0 ");
      }

      file << "\n";
    }

    file << "sub_num\n";

    for (const auto & el : sub_num) {
      if (!el.second.zero()) {
        for (const auto & monomial : el.second.coefs) {
          if (monomial.second != 0) {
            //std::string tmp_str = "";
            file << el.first << " ";
            //tmp_str += std::to_string(el.first) + " ";

            for (uint32_t i = 0; i != n; ++i) {
              file << monomial.first[i] << " ";
            }

            file << " \n";
          }
        }
      }
    }

    file << "sub_den\n";

    for (const auto & el : sub_den) {
      if (!el.second.zero()) {
        for (const auto & monomial : el.second.coefs) {
          if (monomial.second != 0) {
            file << el.first << " ";

            for (uint32_t i = 0; i != n; ++i) {
              file << monomial.first[i] << " ";
            }

            file << " \n";
          }
        }
      }
    }

    file << "zero_degs_num\n";

    if (!zero_degs_num.empty()) {
      for (const auto & el : zero_degs_num) {
        file << el << " ";
      }

      file << "\n";
    }

    file << "zero_degs_den\n";

    if (!zero_degs_den.empty()) {
      for (const auto & el : zero_degs_den) {
        file << el << " ";
      }

      file << "\n";
    }

    file << "g_ni\n";

    for (const auto & el : g_ni) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " " ;
      }

      file << el.second.numerator.get_str() << " " << el.second.denominator.get_str() << "\n";
    }

    file << "g_di\n";

    for (const auto & el : g_di) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " ";
      }

      file << el.second.numerator.get_str() << " " << el.second.denominator.get_str() << "\n";
    }

    file << "combined_ni\n";

    for (const auto & el : combined_ni) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " ";
      }

      file << el.second.get_str() << "\n";
    }

    file << "combined_di\n";

    for (const auto & el : combined_di) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " ";
      }

      file << el.second.get_str() << "\n";
    }

    file << "combined_primes_ni\n";

    for (const auto & el : combined_primes_ni) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " ";
      }

      file << el.second.get_str() << "\n";
    }

    file << "combined_primes_di\n";

    for (const auto & el : combined_primes_di) {
      for (uint32_t i = 0; i != n; ++i) {
        file << el.first[i] << " ";
      }

      file << el.second.get_str() << "\n";
    }

    file << "interpolations\n" << interpolations << "\n";

    file.close();

    std::string file_name_new = "ff_save/states/" + tag + "_" + std::to_string(prime_number) + ".gz";

    if (prime_number > 0) {
      std::string old_file_name = "ff_save/states/" + tag + "_" + std::to_string(prime_number - 1) + ".gz";
      std::rename(old_file_name.c_str(), file_name_new.c_str());
    }

    std::rename(file_name.c_str(), file_name_new.c_str());
  }

  std::pair<bool, uint32_t> RatReconst::start_from_saved_file(const std::string& file_name) {
    std::string line;
    std::ifstream ifile(file_name.c_str());
    igzstream file(file_name.c_str());
    bool first = true;
    {
      std::lock_guard<std::mutex> lock_status(mutex_status);
      prime_number = parse_prime_number(file_name) + 1;
    }
    check_interpolation = true;
    bool tmp_need_shift = false;
    from_save_state = true;

    bool is_zero = false;

    std::vector<uint32_t> individual_degrees_num {};
    std::vector<uint32_t> individual_degrees_den {};

    if (ifile.is_open()) {
      ifile.close();

      while (std::getline(file, line)) {
        if (first) {
          first = false;

          if (line == "ZERO")
            is_zero = true;
          else if (line == "tag_name") {
            std::getline(file, line);
            tag_name = line;
            std::getline(file, line);
            std::getline(file, line);
            normalize_to_den = std::stoi(line);
            std::getline(file, line);
            std::getline(file, line);

            if (line != "individual_degrees_den") {
              individual_degrees_num = parse_vector_32(line);
              std::getline(file, line);
              std::getline(file, line);
              individual_degrees_den = parse_vector_32(line);

              is_set_individual_degree_bounds = true;
              individual_degree_bounds = std::vector<std::pair<uint32_t, uint32_t>> (n);

              for (size_t i = 0; i != n; ++i) {
                individual_degree_bounds[i] = std::make_pair(individual_degrees_num[i], individual_degrees_den[i]);
              }
            }

            file.close();
            std::lock_guard<std::mutex> lock_status(mutex_status);
            prime_number = 0;
            check_interpolation = false;

            return std::make_pair(true, 0);
          } else if (line != "combined_prime") {
            ERROR_MSG("Wrong input format! Has to start with 'combined_prime'!");
            file.close();
            std::exit(EXIT_FAILURE);
          }

          curr_parsed_variable = COMBINED_PRIME;
          parsed_variables[COMBINED_PRIME] = true;
        } else if (is_zero) {
          if (prime_number > 2) {
            std::lock_guard<std::mutex> lock_status(mutex_status);
            new_prime = false;
            done = true;
            g_ni[std::vector<uint32_t>(n)] = RationalNumber(0, 1);
            g_di[std::vector<uint32_t>(n)] = RationalNumber(1, 1);
          } else {
            if (line == "1") {
              set_safe_interpolation();
            }

            check_interpolation = false;
            zero_counter = prime_number;
            {
              std::lock_guard<std::mutex> lock(mutex_status);
              std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
              new_prime = true;
            }
          }

          std::getline(file, line);
          tag_name = line;
          file.close();
          return std::make_pair(false, prime_number);
        } else {
          if (line == "is_done") {
            curr_parsed_variable = IS_DONE;
            parsed_variables[IS_DONE] = true;
          } else if (line == "tag_name") {
            curr_parsed_variable = TAG_NAME;
            parsed_variables[TAG_NAME] = true;
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
          } else if (line == "sub_num") {
            curr_parsed_variable = SUB_NUM;
            parsed_variables[SUB_NUM] = true;
          } else if (line == "sub_den") {
            curr_parsed_variable = SUB_DEN;
            parsed_variables[SUB_DEN] = true;
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
          } else if (line == "combined_primes_ni") {
            curr_parsed_variable = COMBINED_PRIMES_NI;
            parsed_variables[COMBINED_PRIMES_NI] = true;
          } else if (line == "combined_primes_di") {
            curr_parsed_variable = COMBINED_PRIMES_DI;
            parsed_variables[COMBINED_PRIMES_DI] = true;
          } else if (line == "interpolations") {
            curr_parsed_variable = INTERPOLATIONS;
            parsed_variables[INTERPOLATIONS] = true;
          } else if (line == "individual_degrees_num") {
            curr_parsed_variable = INDIVIDUAL_DEGREES_NUM;
          } else if (line == "individual_degrees_den") {
            curr_parsed_variable = INDIVIDUAL_DEGREES_DEN;
          } else {
            switch (curr_parsed_variable) {
              case COMBINED_PRIME: {
                std::lock_guard<std::mutex> lock_status(mutex_status);
                combined_prime = mpz_class(line);
                break;
              }

              case TAG_NAME: {
                tag_name = line;
                break;
              }

              case IS_DONE: {
                std::lock_guard<std::mutex> lock_status(mutex_status);
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

              case INDIVIDUAL_DEGREES_NUM: {
                individual_degrees_num = parse_vector_32(line);
                break;
              }

              case INDIVIDUAL_DEGREES_DEN: {
                individual_degrees_den = parse_vector_32(line);
                break;
              }

              case NEED_PRIME_SHIFT: {
                tmp_need_shift = std::stoi(line);
                std::lock_guard<std::mutex> lock_statics(mutex_statics);

                if (!is_done() && tmp_need_shift && singular_system_set.find(prime_number) == singular_system_set.end())
                  singular_system_set.emplace(prime_number);

                break;
              }

              case NORMALIZER_DEG: {
                normalizer_deg = parse_vector_32(line);
                std::lock_guard<std::mutex> lock_status(mutex_status);

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
                  std::vector<uint32_t> tmp_vec = parse_vector_32(line);
                  std::lock_guard<std::mutex> lock_statics(mutex_statics);

                  for (uint32_t i = 0; i != n; ++i) {
                    if (tmp_vec[i] != 0 && shift[i] == 0)
                      shift[i] = FFInt(get_rand_64());
                    else if (tmp_vec[i] == 0 && shift[i] != 0)
                      shift[i] = 0;
                  }
                }

                break;
              }

              case SUB_NUM: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                uint32_t tmp_deg = parse_vector_32(line, 1)[0];
                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);

                if (sub_num.find(tmp_deg) == sub_num.end())
                  sub_num[tmp_deg] = PolynomialFF(n, {{tmp_vec, 1}});
                else
                  sub_num[tmp_deg] += PolynomialFF(n, {{tmp_vec, 1}});

                break;
              }

              case SUB_DEN: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                uint32_t tmp_deg = parse_vector_32(line, 1)[0];
                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);

                if (sub_den.find(tmp_deg) == sub_den.end())
                  sub_den[tmp_deg] = PolynomialFF(n, {{tmp_vec, 1}});
                else
                  sub_den[tmp_deg] += PolynomialFF(n, {{tmp_vec, 1}});

                break;
              }

              case ZERO_DEGS_NUM: {
                std::vector<uint32_t> tmp_vec = parse_vector_32(line);

                for (const auto & el : tmp_vec) {
                  zero_degs_num.emplace(el);
                }

                break;
              }

              case ZERO_DEGS_DEN: {
                std::vector<uint32_t> tmp_vec = parse_vector_32(line);

                for (const auto & el : tmp_vec) {
                  zero_degs_den.emplace(el);
                }

                break;
              }

              case G_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                g_ni.emplace(std::make_pair(tmp_vec, parse_rational_number(line)));

                break;
              }

              case G_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                g_di.emplace(std::make_pair(tmp_vec, parse_rational_number(line)));

                break;
              }

              case COMBINED_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                combined_ni.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case COMBINED_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                combined_di.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case COMBINED_PRIMES_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                combined_primes_ni.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case COMBINED_PRIMES_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                std::vector<uint32_t> tmp_vec = parse_vector_32(line, n);
                combined_primes_di.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case INTERPOLATIONS: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'normalizer_deg' first.");
                  std::exit(EXIT_FAILURE);
                }

                interpolations = std::stoi(line);

                break;
              }
            }
          }
        }
      }

      size_t counter = 0;
      for (const auto & el : parsed_variables) {
        if (!el) {
          ERROR_MSG("Incomplete input file '" + file_name + "'! It cannot be used to resume a run. Entry '" + save_variables_string[counter] + "' is missing.");
          std::exit(EXIT_FAILURE);
        }
        ++counter;
      }

      file.close();

      if (individual_degrees_num.size() != 0) {
        is_set_individual_degree_bounds = true;
        individual_degree_bounds = std::vector<std::pair<uint32_t, uint32_t>> (n);

        for (size_t i = 0; i != n; ++i) {
          individual_degree_bounds[i] = std::make_pair(individual_degrees_num[i], individual_degrees_den[i]);
        }
      }

      for (const auto & el : combined_ni) add_non_solved_num(el.first);

      for (const auto & el : combined_di) add_non_solved_den(el.first);

      {
        std::lock_guard<std::mutex> lock(mutex_status);
        std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
        new_prime = true;
      }

      if (prime_number >= interpolations) {
        std::lock_guard<std::mutex> lock(mutex_status);
        num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
      }

      for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<FFInt> {};

      for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<FFInt> {};

      return std::make_pair(tmp_need_shift, prime_number);
    } else {
      ERROR_MSG("The file '" + file_name + "' could not be found!");
      std::exit(EXIT_FAILURE);
    }
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

    for (const auto & el : non_solved_degs_num) {
      coef_mat_num[el.first] = std::vector<FFInt> {};
      const uint32_t size = el.second.size();

      if (size > max_num_coef_num.second || (size == max_num_coef_num.second && el.first > max_num_coef_num.first))
        max_num_coef_num = std::make_pair(el.first, size);
    }

    for (const auto & el : non_solved_degs_den) {
      coef_mat_den[el.first] = std::vector<FFInt> {};
      const uint32_t size = el.second.size();

      if (size > max_num_coef_den.second || (size == max_num_coef_den.second && el.first > max_num_coef_den.first))
        max_num_coef_den = std::make_pair(el.first, size);
    }

    PolynomialFF zero_poly(n, {{std::vector<uint32_t>(n), 0}});

    // Initialize subtraction terms with zero
    for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_num); ++i) {
      if (sub_num.find(i) == sub_num.end())
        sub_num[i] = zero_poly;
    }

    for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_den); ++i) {
      if (sub_den.find(i) == sub_den.end())
        sub_den[i] = zero_poly;
    }

    std::vector<uint32_t> zero_vec(n);

    for (auto & el : sub_num) {
      if (el.first != max_num_coef_num.first) {
        if (non_solved_degs_num.find(el.first) != non_solved_degs_num.end()) {
          el.second.remove_zero_coefs();
          std::unordered_set<std::vector<uint32_t>, UintHasher> tmp_set;
          std::vector<std::vector<uint32_t>> tmp_vec = non_solved_degs_num[el.first];

          for (const auto & vec_deg : non_solved_degs_num[el.first]) {
            tmp_set.emplace(vec_deg);
          }

          for (const auto & el2 : el.second.coefs) {
            if (tmp_set.find(el2.first) == tmp_set.end() && el2.first != zero_vec)
              tmp_vec.emplace_back(el2.first);
          }

          if (tmp_vec.size() > 0 && tmp_vec.size() < max_num_coef_num.second) {
            dense_solve_degs_num.emplace(el.first);

            non_solved_degs_num[el.first] = tmp_vec;
          }
        } else {
          if (!(!normalize_to_den && el.first == 0) && !el.second.zero()) {
            el.second.remove_zero_coefs();

            if (el.second.coefs.size() < max_num_coef_num.second) {
              dense_solve_degs_num.emplace(el.first);
              std::vector<std::vector<uint32_t>> tmp_vec;

              if (el.first != 0) {
                for (const auto & el2 : el.second.coefs) {
                  if (el2.first != zero_vec)
                    tmp_vec.emplace_back(el2.first);
                }
              } else {
                for (const auto & el2 : el.second.coefs) {
                  tmp_vec.emplace_back(el2.first);
                }
              }

              non_solved_degs_num[el.first] = tmp_vec;
              coef_mat_num[el.first] = std::vector<FFInt> {};
            } else
              non_solved_degs_num[el.first] = {zero_vec};
          }
        }
      }
    }

    for (auto & el : sub_den) {
      if (el.first != max_num_coef_den.first) {
        if (non_solved_degs_den.find(el.first) != non_solved_degs_den.end()) {
          el.second.remove_zero_coefs();
          std::unordered_set<std::vector<uint32_t>, UintHasher> tmp_set;
          std::vector<std::vector<uint32_t>> tmp_vec = non_solved_degs_den[el.first];

          for (const auto & vec_deg : non_solved_degs_den[el.first]) {
            tmp_set.emplace(vec_deg);
          }

          for (const auto & el2 : el.second.coefs) {
            if (tmp_set.find(el2.first) == tmp_set.end() && el2.first != zero_vec)
              tmp_vec.emplace_back(el2.first);
          }

          if (tmp_vec.size() > 0 && tmp_vec.size() < max_num_coef_den.second) {
            dense_solve_degs_den.emplace(el.first);

            non_solved_degs_den[el.first] = tmp_vec;
          }
        } else {
          if (!(normalize_to_den && el.first == 0) && !el.second.zero()) {
            el.second.remove_zero_coefs();

            if (el.second.coefs.size() < max_num_coef_den.second) {
              dense_solve_degs_den.emplace(el.first);
              std::vector<std::vector<uint32_t>> tmp_vec;

              if (el.first != 0) {
                for (const auto & el2 : el.second.coefs) {
                  if (el2.first != zero_vec)
                    tmp_vec.emplace_back(el2.first);
                }
              } else {
                for (const auto & el2 : el.second.coefs) {
                  tmp_vec.emplace_back(el2.first);
                }
              }

              non_solved_degs_den[el.first] = tmp_vec;
              coef_mat_den[el.first] = std::vector<FFInt> {};
            } else
              non_solved_degs_den[el.first] = {zero_vec};
          }
        }
      }
    }

    sub_num = polff_map();
    sub_den = polff_map();

    // Initialize subtraction terms with zero
    for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_num); ++i) {
      sub_num[i] = zero_poly;
    }

    for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_den); ++i) {
      sub_den[i] = zero_poly;
    }

    {
      std::lock_guard<std::mutex> lock(mutex_status);
      num_eqn = shifted_max_num_eqn;
    }
  }

  void RatReconst::reset(bool change_prime) {
    if (change_prime)
      FFInt::set_new_prime(primes()[0]);

    std::lock_guard<std::mutex> lock(mutex_statics);
    shift = std::vector<FFInt> ();
    singular_system_set = std::unordered_set<uint32_t>();
    global_anchor_points.clear();
    curr_shift = std::vector<uint32_t>();
    PolyReconst::reset();
  }

  void RatReconst::set_safe_interpolation() {
    interpolations = 300;
  }

  bool RatReconst::check_if_done(const FFInt& num, const FFInt& ti) {
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      // TODO clean up
      if (!is_singular_system && (singular_system_set.find(prime_number) != singular_system_set.end()
                                  || shift != std::vector<FFInt> (n)) && prime_number >= interpolations) {
        lock_statics.unlock();
        set_singular_system_vars();
      }
    }

    needed_feed_vec.clear();
    sub_num = polff_map();
    sub_den = polff_map();
    solved_degs_num.clear();
    solved_degs_den.clear();

    if (rec_rat_coef()) {
      bool tmp_done = test_guess(num, ti);
      {
        std::lock_guard<std::mutex> lock(mutex_status);
        done = tmp_done;
      }

      if (done) {
        is_singular_system = false;

        if (tag.size() != 0)
          save_state();

        std::lock_guard<std::mutex> lock(mutex_status);
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

    if (n > 1 && prime_number >= interpolations) {
      PolynomialFF zero_poly(n, {{std::vector<uint32_t>(n), 0}});

      // Initialize subtraction terms with zero
      for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_num); ++i) {
        sub_num[i] = zero_poly;
      }

      for (uint32_t i = 0; i <= static_cast<uint32_t>(max_deg_den); ++i) {
        sub_den[i] = zero_poly;
      }
    }

    if (!use_chinese_remainder) use_chinese_remainder = true;

    {
      std::lock_guard<std::mutex> lock(mutex_status);
      new_prime = false;
    }

    if (!is_singular_system && prime_number >= interpolations) {
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

      for (const auto & el : tmp_num) {
        uint32_t deg = 0;

        for (const auto & tmp_deg : el.first) deg += tmp_deg;

        if (solved_degs_num[deg].zero())
          solved_degs_num[deg] = PolynomialFF(n, {{el.first, el.second}});
        else
          solved_degs_num[deg].coefs.emplace(std::make_pair(el.first, el.second));
      }

      for (const auto & el : tmp_den) {
        uint32_t deg = 0;

        for (const auto & tmp_deg : el.first) deg += tmp_deg;

        if (solved_degs_den[deg].zero())
          solved_degs_den[deg] = PolynomialFF(n, {{el.first, el.second}});
        else
          solved_degs_den[deg].coefs.emplace(std::make_pair(el.first, el.second));
      }
    }

    if (prime_number < interpolations) {
      std::lock_guard<std::mutex> lock(mutex_status);
      zi = 1;
    } else { // Get total amount of needed feeds to interpolate this function over the current prime
      std::map<uint32_t, uint32_t> r_map {};

      // Fill the feed vector
      uint32_t last_number_of_terms = 0;
      uint32_t tmp_max_num_eqn = num_eqn;

      if (is_singular_system) {
        // subtraction term for shifted interpolation
        uint32_t sparse_num = 0;
        uint32_t sparse_den = 0;

        // Check when we can remove functions from the sytem of equations
        for (const auto & el : coef_mat_num) {
          bool got_max_size = false;
          uint32_t size = non_solved_degs_num[el.first].size();

          if (size == max_num_coef_num.second) {
            sparse_num ++;
            got_max_size = true;
          }

          if (got_max_size || dense_solve_degs_num.find(el.first) != dense_solve_degs_num.end()) {
            if (r_map.find(size) != r_map.end())
              r_map[size] ++;
            else
              r_map[size] = 1;
          }
        }

        for (const auto & el : coef_mat_den) {
          bool got_max_size = false;
          uint32_t size = non_solved_degs_den[el.first].size();

          if (size == max_num_coef_den.second) {
            sparse_den ++;
            got_max_size = true;
          }

          if (got_max_size || dense_solve_degs_den.find(el.first) != dense_solve_degs_den.end()) {
            if (r_map.find(size) != r_map.end())
              r_map[size] ++;
            else
              r_map[size] = 1;
          }
        }

        bool num_done = false;
        bool den_done = false;

        for (const auto & el : r_map) {
          uint32_t multiplicity;
          uint32_t tmp_num_eqn = tmp_max_num_eqn;

          if (last_number_of_terms == 0)
            multiplicity = el.first;
          else {
            multiplicity = el.first - last_number_of_terms;

            tmp_num_eqn -= r_map[last_number_of_terms];

            if (num_done) {
              tmp_num_eqn -= (non_solved_degs_num.size() - dense_solve_degs_num.size() - sparse_num); // correct for already removed degrees
              num_done = false;
            }

            if (den_done) {
              tmp_num_eqn -= (non_solved_degs_den.size() - dense_solve_degs_den.size() - sparse_den); // correct for already removed degrees
              den_done = false;
            }
          }

          if (!num_done && el.first == max_num_coef_num.second)
            num_done = true;

          if (!den_done && el.first == max_num_coef_den.second)
            den_done = true;

          needed_feed_vec.emplace_back(std::make_pair(multiplicity, tmp_num_eqn));
          last_number_of_terms = el.first;
          tmp_max_num_eqn = tmp_num_eqn;
        }
      } else {
        if (n > 1) {
          // Check when we can remove functions from the sytem of equations
          for (const auto & el : coef_mat_num) {
            uint32_t size = non_solved_degs_num[el.first].size();

            if (r_map.find(size) != r_map.end())
              r_map[size] ++;
            else
              r_map[size] = 1;
          }

          for (const auto & el : coef_mat_den) {
            uint32_t size = non_solved_degs_den[el.first].size();

            if (r_map.find(size) != r_map.end())
              r_map[size] ++;
            else
              r_map[size] = 1;
          }

          for (const auto & el : r_map) {
            uint32_t multiplicity;
            uint32_t tmp_num_eqn = tmp_max_num_eqn;

            if (last_number_of_terms == 0)
              multiplicity = el.first;
            else {
              multiplicity = el.first - last_number_of_terms;

              tmp_num_eqn -= r_map[last_number_of_terms];
            }

            needed_feed_vec.emplace_back(std::make_pair(multiplicity, tmp_num_eqn));
            last_number_of_terms = el.first;
            tmp_max_num_eqn = tmp_num_eqn;
          }
        } else {
          needed_feed_vec.emplace_back(std::make_pair(1, num_eqn));
        }
      }
    }

    // Rewrite needed_feed_vec in terms of zi_orders. Note that only one entry is necessary here
    std::vector<std::pair<uint32_t, uint32_t>> needed_feed_vec_tmp {};
    size_t counter = 1;

    for (const auto& el : needed_feed_vec) {
      uint32_t tmp_mult = el.first;
      for (uint32_t i = 0; i != tmp_mult; ++i) {
        needed_feed_vec_tmp.emplace_back(std::make_pair(counter, el.second));
        ++counter;
      }
    }

    needed_feed_vec_tmp.shrink_to_fit();
    needed_feed_vec.swap(needed_feed_vec_tmp);

    if (from_save_state) {
      from_save_state = false;

      if (prime_number < interpolations) {
        if (parsed_probes.size() > 0) {

          for (const auto & el : parsed_probes) {
            for (const auto & el2 : el.second) {
              queue.emplace(std::make_tuple(el2.first, el2.second, el.first));
            }
          }
        }
      } else {
        // Calculate the difference of already calculated probes and remaining ones
        if (parsed_probes.size() != 0) {
          size_t size = needed_feed_vec.size();

          for (const auto & el : parsed_probes) {
            if (size == 0) {
              break;
            }

            size_t tmp_size = el.second.size();

            if (n > 1) {
              if (el.first[0] <= size) {
                if (tmp_size >= needed_feed_vec[el.first[0] - 1].second) {
                  needed_feed_vec[el.first[0] - 1].second = 0;
                } else {
                  needed_feed_vec[el.first[0] - 1].second -= tmp_size;
                }
              }
            } else {
              if (tmp_size >= needed_feed_vec[0].second) {
                needed_feed_vec.clear();
              } else {
                needed_feed_vec[0].second -= tmp_size;
              }
            }

            for (const auto & el2 : el.second) {
              queue.emplace(std::make_tuple(el2.first, el2.second, el.first));
            }
          }
        }
      }

      parsed_probes.clear();
    }

    return false;
  }

  std::vector<std::pair<uint32_t, uint32_t>> RatReconst::get_needed_feed_vec() {
    std::vector<std::pair<uint32_t, uint32_t>> needed_feed_vec_tmp = std::move(needed_feed_vec);
    needed_feed_vec.clear();

    return needed_feed_vec_tmp;
  }

  std::vector<FFInt> RatReconst::get_anchor_points() {
    return private_anchor_points;
  }

  std::pair<uint32_t, uint32_t> RatReconst::get_max_deg() {
    if (all_shift_max_degs.size() != 0) {
      return std::make_pair(all_shift_max_degs[0], all_shift_max_degs[1]);
    } else if (max_deg_num != -1 && max_deg_den != -1)
      return std::make_pair(max_deg_num, max_deg_den);
    else {
      if (!is_calc_factors)
        WARNING_MSG("Maximum degrees are not known yet.");

      return std::make_pair(0, 0);
    }
  }

  std::string RatReconst::get_tag() {
    return tag;
  }

  std::string RatReconst::get_tag_name() {
    return tag_name;
  }

  void RatReconst::write_food_to_file() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (!is_writing_probes && saved_food.size() != 0) {
      is_writing_probes = true;
      start = std::chrono::high_resolution_clock::now();
      std::queue<std::tuple<FFInt, FFInt, std::vector<uint32_t>>> tmp;
      tmp.swap(saved_food);
      ogzstream file;
      std::string file_name = "ff_save/probes/" + tag + "_" + std::to_string(prime_number) + ".gz";
      lock.unlock();

      file.open(file_name.c_str(), std::ios_base::app);

      while (!tmp.empty()) {
        auto tmp_el = tmp.front();
        tmp.pop();

        for (const auto & el2 : std::get<2>(tmp_el)) {
          file << el2 << " ";
        }

        file << std::get<0>(tmp_el) << " " << std::get<1>(tmp_el) << " \n";
      }

      file.close();
      lock.lock();
      is_writing_probes = false;
    }
  }

  std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> RatReconst::read_in_probes(const std::string& file_name) {
    std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> tmp {};
    std::string line;
    igzstream file(file_name.c_str());
    auto prime_it = parse_prime_number(file_name);

    if (!done && prime_it == prime_number) {
      while (std::getline(file, line)) {
        std::vector<uint32_t> tmp_zi_ord;

        if (n > 1)
          tmp_zi_ord = parse_vector_32(line, n - 1);

        std::vector<FFInt> tmp_probe = parse_vector_FFInt(line, 2);

        tmp[tmp_zi_ord].emplace(tmp_probe[0].n);

        if (prime_it == 0)
          queue.emplace(std::make_tuple(tmp_probe[0], tmp_probe[1], tmp_zi_ord));
        else
          parsed_probes[tmp_zi_ord].emplace_back(std::make_pair(tmp_probe[0].n, tmp_probe[1].n));
      }
    }

    file.close();

    return tmp;
  }

  void RatReconst::set_prime_to_max() {
    std::lock_guard<std::mutex> lock(mutex_status);
    prime_number = 299;
  }

  std::pair<std::vector<std::pair<std::vector<uint32_t>, uint32_t>>, uint32_t> RatReconst::get_zi_orders() const {
    std::lock_guard<std::mutex> lock(mutex_status);
    uint32_t sub_term = 0;
    if (saved_ti.find(curr_zi_order) != saved_ti.end())
      sub_term = saved_ti.at(curr_zi_order).size();
    else if(coef_mat.size() != 0)
      sub_term = coef_mat.size();

    if (zi < 3) {
      return std::make_pair(std::vector<std::pair<std::vector<uint32_t>, uint32_t>> (1, std::make_pair(curr_zi_order, num_eqn - sub_term)), num_eqn);
    } else {
      std::vector<std::pair<std::vector<uint32_t>, uint32_t>> tmp_zi_orders (vandermonde_size, std::make_pair(curr_zi_order, num_eqn - sub_term));
      for (size_t i = 1; i != vandermonde_size; ++i) {
        for (uint32_t tmp_zi = 2; tmp_zi != zi; ++tmp_zi) {
          tmp_zi_orders[i].first[tmp_zi - 2] += i;
        }

        uint32_t tmp_sub_term = 0;

        if (saved_ti.find(tmp_zi_orders[i].first) != saved_ti.end())
          tmp_sub_term = saved_ti.at(tmp_zi_orders[i].first).size();

        tmp_zi_orders[i].second = num_eqn - tmp_sub_term;
      }

      return std::make_pair(tmp_zi_orders, num_eqn);
    }
  }

  void RatReconst::set_individual_degree_bounds(const std::vector<std::pair<uint32_t, uint32_t>>& individual_degree_bounds_) {
    is_set_individual_degree_bounds = true;
    individual_degree_bounds = individual_degree_bounds_;
  }

  std::vector<uint32_t> RatReconst::get_individual_degree_bounds(bool is_den) {
    std::vector<uint32_t> res (n - 1);

    for (size_t i = 1; i != n; ++i) {
      if (is_den) {
        res[i - 1] = individual_degree_bounds[i].second;
      } else {
        res[i - 1] = individual_degree_bounds[i].first;
      }
    }

    return res;
  }
}
