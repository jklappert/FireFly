#include "Logger.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sys/stat.h>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::need_prime_shift = false;
  bool RatReconst::set_singular_system = false;
  ff_pair_map RatReconst::rand_zi;
  std::mutex RatReconst::mutex_statics;

  RatReconst::RatReconst(uint32_t n_) {
    n = n_;
    type = RAT;
    std::unique_lock<std::mutex> lock_status(mutex_status);

    combined_prime = FFInt::p;

    //std::srand(std::time(nullptr));
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (shift.empty()) {
        shift = std::vector<FFInt> (n);

        if (n > 1) {
          for (auto & el : shift) el = get_rand();

          curr_zi_order_num = std::vector<uint32_t> (n - 1, 1);
          curr_zi_order_den = std::vector<uint32_t> (n - 1, 1);
        }
      }
    }

    if (n > 1) {
      curr_zi_order = std::vector<uint32_t> (n - 1, 1);
      lock_status.unlock();
      // add a zero to both polynomials to do arithmetic
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

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t& fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (!done && fed_prime == prime_number)
      queue.emplace_back(std::make_tuple(new_ti, num, feed_zi_ord));
  }

  void RatReconst::interpolate() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (is_interpolating || queue.empty()) return;
    else {
      is_interpolating = true;

      while (!queue.empty()) {
        auto food = queue.front();
        queue.pop_front();
        lock.unlock();
        interpolate(std::get<0>(food), std::get<1>(food), std::get<2>(food));

        while (saved_ti.find(curr_zi_order) != saved_ti.end()) {
          /*
          * If not finished, check if we can use some saved runs
          */
          if (saved_ti.find(curr_zi_order) != saved_ti.end()) {
            std::pair<FFInt, FFInt> key_val = saved_ti.at(curr_zi_order).back();
            saved_ti.at(curr_zi_order).pop_back();
            interpolate(key_val.first, key_val.second, curr_zi_order);
          }
        }

        lock.lock();
      }
    }

    is_interpolating = false;
  }

  void RatReconst::interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint32_t>& feed_zi_ord) {
    if (!done) {

      // Compare if the food is the expected food; if not, store it for later use
      if (feed_zi_ord == curr_zi_order) {
        {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);

          if (!is_singular_system && set_singular_system)
            set_singular_system_vars();
        }

        // first check if we are done. If not start the reconstruction again using
        // the chinese remainder theorem in combining the previous results
        if (new_prime) {
          ti.emplace_back(new_ti);

          if (rec_rat_coef()) {
            bool tmp_done = test_guess(num);
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              done = tmp_done;
            }

            if (done) {
              std::unique_lock<std::mutex> lock(mutex_status);
              combined_di = mpz_map();
              combined_ni = mpz_map();
              combined_prime = 0;
              num_eqn = 0;
              new_prime = false;
              solved_den = PolynomialFF();
              solved_num = PolynomialFF();
              coef_mat_num = std::unordered_map<uint32_t, std::vector<std::pair<FFInt, uint32_t>>> ();
              coef_mat_den = std::unordered_map<uint32_t, std::vector<std::pair<FFInt, uint32_t>>> ();
              curr_zi_order.clear();
              use_chinese_remainder = false;
              return;
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
          ti.pop_back();

          if (!is_singular_system) {
            ff_map tmp_num {};
            tmp_num.reserve(g_ni.size() + 1);
            ff_map tmp_den {};
            tmp_den.reserve(g_di.size() + 1);
            tmp_num.emplace(std::make_pair(std::vector<uint32_t> (n), 0));
            tmp_den.emplace(std::make_pair(std::vector<uint32_t> (n), 0));

            for (const auto & el : g_ni) {
              tmp_num.emplace(std::make_pair(el.first, FFInt(el.second.numerator) / FFInt(el.second.denominator)));
            }

            for (const auto & el : g_di) {
              tmp_den.emplace(std::make_pair(el.first, FFInt(el.second.numerator) / FFInt(el.second.denominator)));
            }

            solved_num = PolynomialFF(n, tmp_num);
            solved_den = PolynomialFF(n, tmp_den);
          }
        }

        // basic reconstruction algorithm, check if reconstructed function is equal
        // to numeric input and calculate coefficients a_i, check chinese chinese remainder
        // theorem
        {
          std::unique_lock<std::mutex> lock(mutex_status);

          if (prime_number == 0) zi = 1;
        }

        if (max_deg_num == -1) {
          ti.emplace_back(new_ti);
          const uint32_t i = ti.size() - 1;

          if (i == 0) {
            ai.emplace_back(num);
          } else {
            if (num == comp_fyi(i - 1, i - 1, ti.back())) check = true;

            if (!check)
              ai.emplace_back(comp_ai(i, i, num));
          }
        } else {
          if (coef_mat.empty())
            coef_mat.reserve(num_eqn);

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

            if (n > 1) {
              yis = get_rand_zi_vec(curr_zi_order);
            }

            yis.emplace(yis.begin(), 1);

            // build Gauss system for univariate reconstruction needed for
            // multivariate rational functions
            if (prime_number == 0)
              build_uni_gauss(tmp_ti, tmp_num, yis);
            else
              build_homogenized_multi_gauss(tmp_ti, tmp_num, yis);

            if (coef_mat.size() == num_eqn) {
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

            ti.pop_back();

            canonical = construct_canonical();
            PolynomialFF numerator = PolynomialFF(1, canonical.first);
            PolynomialFF denominator = PolynomialFF(1, canonical.second);

            //TODO catch new shift
            if (n > 1 && denominator.min_deg()[0] > 0) {
              INFO_MSG("No constant term in denominator! Trying again with new paramter shift...");

              for (uint32_t j = 0; j < n; j++) {
                shift[j] = get_rand();
              }

              {
                std::unique_lock<std::mutex> lock(mutex_status);
                done = false;
              }

              ai.clear();
              ti.clear();
              return;
            }

            max_deg_num = numerator.max_deg()[0];
            max_deg_den = denominator.max_deg()[0];

            if (n == 1)
              min_deg_den_vec = denominator.min_deg();

            curr_deg_num = max_deg_num;

            if (max_deg_den > 0)
              curr_deg_den = max_deg_den;

            FFInt equializer = FFInt(1) / denominator.coefs[denominator.min_deg()];

            canonical.first = (numerator * equializer).coefs;
            canonical.second = (denominator * equializer).coefs;

            // set number of equations needed for univariate rational function
            // reconstruction needed for multivariate polynomial feed
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              num_eqn = max_deg_den + max_deg_num + 1
                        - tmp_solved_coefs_num - tmp_solved_coefs_den;
            }
            ai.clear();
            ti.clear();
          } else if (prime_number == 0)
            canonical = solve_gauss();
          else
            canonical = solve_homogenized_multi_gauss();

          if (n == 1) {
            std::pair<mpz_map, mpz_map> tmp;
            tmp.first = convert_to_mpz(canonical.first);
            tmp.second = convert_to_mpz(canonical.second);
            combine_primes(tmp);
            saved_ti.clear();
            std::unique_lock<std::mutex> lock(mutex_status);
            prime_number++;
            queue.clear();
            new_prime = true;
            return;
          } else if (prime_number == 0) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi = curr_zi;
            }

            ff_map num_coef = canonical.first;
            ff_map den_coef = canonical.second;
            std::vector<uint32_t> zero_deg(n);
            ff_map zero_mon = {{zero_deg, 0}};

            // save the current results to the map to access them later
            for (int i = 0; i <= (int)(max_deg_num - tmp_solved_coefs_num); i++) {
              if (first_run) {
                {
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  PolyReconst rec(n - 1, (uint32_t) i, true);

                  if (rec.is_rand_zi_empty()) {
                    std::vector<FFInt> anchor_points {};

                    for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
                      anchor_points.emplace_back(rand_zi[std::make_pair(tmp_zi, 1)]);
                    }

                    rec.set_anchor_points(anchor_points);
                  }

                  coef_n.emplace(std::make_pair((uint32_t) i, std::move(rec)));
                }
              }

              if (i <= curr_deg_num) {
                // this saves some memory since we only need one numerical value
                // for the constant coefficient
                if (i == 0 && first_run) {
                  std::vector<uint32_t> key = {(uint32_t) i, zi};
                  uint32_t sub_count = sub_num[i].size();
                  saved_num_num[curr_zi_order][key] = std::make_pair(num_coef[ {(uint32_t) i}], sub_count);
                } else {
                  if (curr_zi_order[zi - 2] < (uint32_t) i + 3) {
                    uint32_t sub_count = sub_num[i].size();
                    std::vector<uint32_t> key = {(uint32_t) i, zi};
                    saved_num_num[curr_zi_order][key] = std::make_pair(num_coef[ {(uint32_t) i}], sub_count);
                  }
                }
              }
            }

            for (uint32_t i = 1; i <= max_deg_den - tmp_solved_coefs_den; i++) {
              if (first_run) {
                {
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  PolyReconst rec(n - 1, i, true);

                  if (rec.is_rand_zi_empty()) {
                    std::vector<FFInt> anchor_points {};

                    for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi ++)
                      anchor_points.emplace_back(rand_zi[std::make_pair(tmp_zi, 1)]);

                    rec.set_anchor_points(anchor_points);
                  }

                  coef_d.emplace(std::make_pair(i, std::move(rec)));
                }
              }

              if ((int) i <= curr_deg_den) {
                if (curr_zi_order[zi - 2] < (uint32_t) i + 3) {
                  uint32_t sub_count = sub_den[i].size();
                  std::vector<uint32_t> key = {i, zi};
                  saved_num_den[curr_zi_order][key] = std::make_pair(den_coef[ {i}], sub_count);
                }
              }
            }

            if (first_run) first_run = false;

            uint32_t zi_num = 0;

            if (curr_deg_num > 0) zi_num = coef_n[curr_deg_num].get_zi() + 1;

            uint32_t zi_den = 0;

            if (curr_deg_den > 0) zi_den = coef_d[curr_deg_den].get_zi() + 1;

            // reconstruct the numerator
            if (curr_deg_num >= 0) {
              PolyReconst rec = coef_n[curr_deg_num];

              if (rec.get_zi_order() == std::vector<uint32_t>(curr_zi_order.begin(), curr_zi_order.end())) {
                auto res = feed_poly(curr_deg_num, max_deg_num, coef_n, rec,
                                     saved_num_num, sub_num, true);
                curr_deg_num = std::get<0>(res);
                zi_num = std::get<1>(res);
                curr_zi_order_num = std::get<2>(res);
              }
            }

            // reconstruct the denominator
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

              if (rec.get_zi_order() == std::vector<uint32_t>(curr_zi_order.begin(), curr_zi_order.end())) {
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

                if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  numerator += res;
              }

              for (auto & el : coef_d) {
                PolynomialFF res = el.second.get_result_ff();

                if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  denominator += res;
              }

              FFInt terminator = 0;
              // check if one can normalize to a single univariate degree, if so
              // set the corresponding degree in equalizer_degree. Check constants
              // first and proceed with denominator/numerator.
              // if the denominator is just a constant, there is no corresponding
              // PolyReconst object. Thus we set the minimal degree to a zero tuple

              if (const_den != 1) {
                ff_map dummy_map {};
                terminator = FFInt(1) - const_den;
                dummy_map.emplace(std::make_pair(std::vector<uint32_t> (n, 0), terminator));
                denominator = denominator + PolynomialFF(n, dummy_map);
              } else if (numerator.coefs.find(std::vector<uint32_t> (n, 0)) != numerator.coefs.end()) {
                terminator = numerator.coefs[std::vector<uint32_t> (n, 0)];
              } else {
                for (const auto & el : denominator.coefs) {
                  add_non_solved_den(el.first);

                  if (min_deg_den_vec.empty())
                    min_deg_den_vec = el.first;
                  else if(a_grt_b(min_deg_den_vec, el.first))
                    min_deg_den_vec = el.first;
                }

                for (const auto & candidate : non_solved_degs_den) {
                  if (candidate.second.size() == 1) {
                    terminator = denominator.coefs[candidate.second[0]];
                    break;
                  }
                }

                if (terminator.n == 0) {
                  for (const auto & el : numerator.coefs) {
                    add_non_solved_num(el.first);
                  }


                  for (const auto & candidate : non_solved_degs_num) {
                    if (candidate.second.size() == 1) {
                      terminator = numerator.coefs[candidate.second[0]];
                      break;
                    }
                  }
                }

                if (terminator.n == 0) {
                  // normalize to the minimal degree of the denominator
                  terminator = denominator.coefs[min_deg_den_vec];

                  is_singular_system = true;
                }
              }

              coef_n.clear();
              coef_d.clear();
              curr_zi_order_num.clear();
              curr_zi_order_den.clear();

              // normalize
              FFInt equializer = FFInt(1) / terminator;

              numerator = numerator * equializer;
              denominator = denominator * equializer;

              std::pair<mpz_map, mpz_map> tmp;
              tmp.first = convert_to_mpz(numerator.coefs);
              tmp.second = convert_to_mpz(denominator.coefs);

              combine_primes(tmp);

              std::unique_lock<std::mutex> lock(mutex_status);
              prime_number++;
              queue.clear();
              saved_ti.clear();
              std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
              curr_zi = 2;
              zi = 2;
              new_prime = true;
            } else if (zi > 0) {
              // set new random
              for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
                auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 2]);
                std::unique_lock<std::mutex> lock_statics(mutex_statics);
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
                  solved_num += solve_transposed_vandermonde(non_solved_degs_num[key], coef_mat_num[key]);
                  non_solved_degs_num.erase(key);
                  coef_mat_num.erase(key);
                }
              }

              for (const auto & sol : canonical.second) {
                uint32_t key = sol.first[0];

                if (coef_mat_den[key].size() < non_solved_degs_den[key].size())
                  coef_mat_den[key].emplace_back(std::make_pair(sol.second, 0));

                if (coef_mat_den[key].size() == non_solved_degs_den[key].size()) {
                  solved_den += solve_transposed_vandermonde(non_solved_degs_den[key], coef_mat_den[key]);
                  non_solved_degs_den.erase(key);
                  coef_mat_den.erase(key);
                }
              }
            } else {
              uint32_t tmp_deg = curr_deg_num;

              if (curr_deg_num > -1) {
                for (uint32_t i = 0; i <= tmp_deg; i++) {
                  if (coef_mat_num.find(i) != coef_mat_num.end()) {
                    if (coef_mat_num[i].size() < non_solved_degs_num[i].size()) {
                      uint32_t sub_count = sub_num[i].size();
                      coef_mat_num[i].emplace_back(std::make_pair(canonical.first[ {i}], sub_count));
                    }

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

              if (curr_deg_den > 0) {
                tmp_deg = curr_deg_den;

                for (uint32_t i = 1; i <= tmp_deg; i++) {
                  if (coef_mat_den.find(i) != coef_mat_den.end()) {
                    if (coef_mat_den[i].size() < non_solved_degs_den[i].size()) {
                      uint32_t sub_count = sub_den[i].size();
                      coef_mat_den[i].emplace_back(std::make_pair(canonical.second[ {i}], sub_count));
                    }

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
                FFInt terminator = solved_den.coefs[min_deg_den_vec];
                FFInt equializer = FFInt(1) / terminator;

                for (const auto & el : g_ni) {
                  solved_num.coefs.erase(el.first);
                }

                for (const auto & el : g_di) {
                  solved_den.coefs.erase(el.first);
                }

                solved_num = solved_num * equializer;
                solved_den = solved_den * equializer;
              }

              solved_num.coefs.erase(std::vector<uint32_t> (n, 0));
              solved_den.coefs.erase(std::vector<uint32_t> (n, 0));

              std::pair<mpz_map, mpz_map> tmp;
              tmp.first = convert_to_mpz(solved_num.coefs);
              tmp.second = convert_to_mpz(solved_den.coefs);
              combine_primes(tmp);
              {
                std::unique_lock<std::mutex> lock(mutex_status);
                prime_number++;
                queue.clear();
                saved_ti.clear();
                std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
                new_prime = true;
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
                std::transform(curr_zi_order.begin(), curr_zi_order.end(),
                curr_zi_order.begin(), [](uint32_t x) {return x + 1;});

                for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
                  auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 2]);
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  rand_zi.emplace(std::make_pair(key, rand_zi[std::make_pair(tmp_zi, 1)].pow(key.second)));
                }

                if (!is_singular_system)
                  num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
                else {
                  num_eqn = max_deg_den + max_deg_num + 1
                            - tmp_solved_coefs_num - tmp_solved_coefs_den;
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
          while (!rec.is_new_prime()) {
            FFInt sub_num = 0;

            if (curr_deg != (int)max_deg && sub_count < sub_save[curr_deg].size()) {
              yis.emplace(yis.begin(), 1);
              sub_num = sub_save[curr_deg][sub_count].calc(yis);
              yis.erase(yis.begin());
            }

            rec.feed(yis, food - sub_num);
          }
        } else {
          FFInt sub_num = 0;

          if (curr_deg != (int)max_deg && sub_count < sub_save[curr_deg].size()) {
            yis.emplace(yis.begin(), 1);
            sub_num = sub_save[curr_deg][sub_count].calc(yis);
            yis.erase(yis.begin());
          }

          rec.feed(yis, food - sub_num);
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
          std::vector<uint32_t> zero_deg(n);
          PolynomialFF zero_poly(n, {{zero_deg, 0}});

          for (int i = 0; i < curr_deg; i++) {
            if (sub_save[(uint32_t)i].size() == 0)
              sub_save[(uint32_t)i] = {zero_poly};
            else
              sub_save[(uint32_t)i].emplace_back(zero_poly);
          }

          PolynomialFF res = rec.get_result_ff();

          // check if the polynomial is zero which we then can omit for further
          // calculations
          if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0)) {
            std::vector<FFInt> tmp_shift;
            {
              std::unique_lock<std::mutex> lock_statics(mutex_statics);
              tmp_shift = shift;
            }
            PolynomialFF sub_pol = rec.get_result_ff().add_shift(tmp_shift);
            sub_pol -= rec.get_result_ff();

            for (auto & el : sub_pol.coefs) {
              int tmp_deg = 0;

              for (const auto & deg : el.first) tmp_deg += deg;

              if (tmp_deg < curr_deg) {
                for (auto & tmp_sub : sub_save[(uint32_t)tmp_deg]) {
                  tmp_sub += PolynomialFF(n, {{el.first, el.second}});
                }
              }
            }
          }

          if (!is_num) {
            std::vector<FFInt> tmp_yis(n, 0);
            const_den += sub_save[0].back().calc(tmp_yis);
          }
        }

        /*
         * Remove already solved coefficients from Gauss eliminiation
         */
        curr_deg--;

        if (is_num)
          tmp_solved_coefs_num ++;
        else {
          if (curr_deg > -1)
            tmp_solved_coefs_den ++;

          if (curr_deg == 0)
            curr_deg = -1;
        }

        {
          std::unique_lock<std::mutex> lock(mutex_status);
          num_eqn = max_deg_den + max_deg_num + 1
                    - tmp_solved_coefs_num - tmp_solved_coefs_den;
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

  void RatReconst::combine_primes(std::pair<mpz_map, mpz_map>& tmp) {
    std::vector<uint32_t> tmp_deg_num {};
    std::vector<uint32_t> tmp_deg_den {};

    if (is_singular_system) {
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

    if (!use_chinese_remainder) {
      saved_num_num = ff_map_map();
      saved_num_den = ff_map_map();
      coef_d = std::unordered_map<uint32_t, PolyReconst>();
      coef_n = std::unordered_map<uint32_t, PolyReconst>();
      combined_ni = tmp.first;
      combined_di = tmp.second;

      // if the coefficient is not a rational number thus divided by 1,
      // it will not change in the next run and can be omitted to save
      // numerical runs
      mpz_map combined_ni_back = combined_ni;

      for (auto & c_ni : combined_ni_back) {
        RationalNumber rn_wang;
        bool wang;

        auto res = get_rational_coef(c_ni.second, combined_prime);
        wang = res.first;

        if (res.first)
          rn_wang = res.second;

        if (wang && rn_wang.numerator == c_ni.second && rn_wang.denominator == 1) {
          remove_ni(c_ni.first, rn_wang);
          continue;
        }

        add_non_solved_num(c_ni.first);
      }

      mpz_map combined_di_back = combined_di;

      for (auto & c_di : combined_di_back) {
        RationalNumber rn_wang;
        bool wang;

        auto res = get_rational_coef(c_di.second, combined_prime);
        wang = res.first;

        if (wang)
          rn_wang = res.second;


        if (wang && rn_wang.numerator == c_di.second && rn_wang.denominator == 1) {
          remove_di(c_di.first, rn_wang);
          continue;
        }

        add_non_solved_den(c_di.first);
      }

      if (is_singular_system) {
        check_for_solved_degs(tmp_deg_num, true);

        if (is_singular_system)
          check_for_solved_degs(tmp_deg_den, false);
      }

      combined_ni_back.clear();
      combined_di_back.clear();
      tmp_solved_coefs_num = 0;
      tmp_solved_coefs_den = 0;
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

      combined_ni_back.clear();
      combined_di_back.clear();
      combined_prime_back = 0;
    }

    const_den = 0;

    if (is_singular_system) {
      tmp_solved_coefs_den = 0;
      tmp_solved_coefs_num = 0;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        need_prime_shift = true;
      }

      for (const auto & el : g_ni) add_non_solved_num(el.first);

      for (const auto & el : g_di) add_non_solved_den(el.first);

      curr_deg_num = max_deg_num;
      curr_deg_den = max_deg_den;

      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 1;
    } else {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
    }

    for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

    for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

    // Check if the state should be written out after this prime
    if (tag.size() > 0) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      std::ofstream file;
      std::string file_name = std::string("ff_save/") + tag + std::string("_") + std::to_string(prime_number) + std::string(".txt");
      file.open(file_name.c_str());
      file << "combined_prime\n" << combined_prime.get_str() << "\n";
      file << "max_deg_num\n" << max_deg_num << "\n";
      file << "max_deg_den\n" << max_deg_den << "\n";
      file << "need_prime_shift\n" << need_prime_shift << "\n";
      file << "min_deg_den_vec\n";
      std::string tmp_vec = "";

      for (const auto & deg : min_deg_den_vec) {
        tmp_vec += std::to_string(deg) + std::string(" ");
      }

      tmp_vec.substr(0, tmp_vec.size() - 1);
      tmp_vec += std::string("\n");
      file << tmp_vec;
      file << "g_ni\n";

      for (const auto & el : g_ni) {
        std::string tmp_entry = "";

        for (const auto & deg : el.first) tmp_entry += std::to_string(deg) + std::string(" ");

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

      file.close();

      if (prime_number > 0) {
        std::string old_file_name = std::string("ff_save/") + tag + std::string("_") + std::to_string(prime_number - 1) + std::string(".txt");

        if (std::remove(old_file_name.c_str()) != 0)
          WARNING_MSG("The previously saved file could not be deleted.");
      }
    }
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
        numerator.sort();
        denominator.sort();

        result = RationalFunction(numerator, denominator);

        RationalNumber first_coef = result.denominator.coefs[0].coef;

        if (first_coef.numerator != 1 || first_coef.denominator != 1) result = normalize(result);
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

  FFInt RatReconst::comp_ai(int i, int ip, const FFInt& num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt ai_i = comp_ai(i, ip - 1, num);
      return (ti[i] - ti[ip - 1]) / (ai_i - ai[ip - 1]);
    }
  }

  FFInt RatReconst::comp_fyi(uint32_t i, uint32_t ip, const FFInt& y) {
    if (ip == 0) {
      return ai[i];
    } else {
      return ai[i - ip] + (-ti[i - ip] + y) / comp_fyi(i, ip - 1, y);
    }
  }

  std::pair<ff_map, ff_map> RatReconst::solve_gauss() {
    std::vector<FFInt> results = solve_gauss_system(num_eqn, coef_mat);
    coef_mat.clear();

    // Bring result in canonical form
    ff_map numerator;
    ff_map denominator;

    const std::vector<uint32_t> min_power = {0};
    denominator.emplace(std::make_pair(std::move(min_power), FFInt(1)));

    int terms_num = max_deg_num - tmp_solved_coefs_num;

    if (terms_num == -1) {
      numerator.emplace(std::make_pair(std::move(min_power), FFInt(1)));
    } else {
      for (int i = 0; i <= terms_num; i ++) {
        std::vector<uint32_t> power = { (uint32_t) i};
        numerator.emplace(std::make_pair(std::move(power), results[i]));
      }
    }

    for (uint32_t i = 1; i <= max_deg_den - tmp_solved_coefs_den; i ++) {
      std::vector<uint32_t> power = {i};
      denominator.emplace(std::make_pair(std::move(power), results[i + terms_num]));
    }

    return std::make_pair(numerator, denominator);
  }

  std::pair<ff_map, ff_map> RatReconst::solve_homogenized_multi_gauss() {
    std::vector<FFInt> results = solve_gauss_system(num_eqn, coef_mat);
    coef_mat.clear();

    // Bring result in canonical form
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
      int terms_num = max_deg_num - tmp_solved_coefs_num;

      if (terms_num > 0) {
        for (int i = 0; i <= terms_num; i++) {
          if (non_solved_degs_num.find((uint32_t)i) != non_solved_degs_num.end()) {
            std::vector<uint32_t> power = {(uint32_t)i};
            numerator.emplace(std::make_pair(std::move(power), results[counter]));
          }

          counter ++;
        }
      }

      for (uint32_t i = 1; i <= max_deg_den - tmp_solved_coefs_den; i++) {
        if (non_solved_degs_den.find(i) != non_solved_degs_den.end()) {
          std::vector<uint32_t> power = {i};
          denominator.emplace(std::make_pair(std::move(power), results[counter]));
        }

        counter ++;
      }
    }

    return std::make_pair(numerator, denominator);
  }

  std::pair<ff_map, ff_map> RatReconst::construct_canonical() {
    if (ai.size() == 1) {
      ff_map numerator_ff;
      std::vector<uint32_t> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[0]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(numerator_ff, denominator_ff);
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(1);
      FFInt mti = -ti[0];
      return std::make_pair((r.first * ai[0] + r.second * mti + r.second.mul(1)).coefs,
                            r.first.coefs);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::iterate_canonical(uint32_t i) {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(i + 1);
      FFInt mti = -ti[i];
      return std::pair<PolynomialFF, PolynomialFF> (fnp1.first * ai[i] + fnp1.second.mul(1) + fnp1.second * mti,
                                                    fnp1.first);
    } else {
      ff_map numerator_ff;
      std::vector<uint32_t> zero_deg = {0};
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

  bool RatReconst::test_guess(const FFInt& num) {
    ff_map g_ff_ni = convert_to_ffint(g_ni);
    ff_map g_ff_di = convert_to_ffint(g_di);
    PolynomialFF g_ny(n, g_ff_ni);
    PolynomialFF g_dy(n, g_ff_di);
    std::vector<FFInt> yis = std::vector<FFInt> (n);
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yis[0] = ti[0] + shift[0];

      for (uint32_t i = 1; i < n; i++) {
        yis[i] = ti[0] * rand_zi[std::make_pair(i + 1, curr_zi_order[i - 1])] + shift[i];
      }
    }

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;
  }

  void RatReconst::remove_ni(const std::vector<uint32_t>& deg_vec, RationalNumber& rn) {
    g_ni[deg_vec] =  rn;
    combined_ni.erase(deg_vec);
  }

  void RatReconst::remove_di(const std::vector<uint32_t>& deg_vec, RationalNumber& rn) {
    g_di[deg_vec] =  rn;
    combined_di.erase(deg_vec);
  }

  RatReconst::RatReconst(const RatReconst& other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    first_run = other.first_run;
    coef_mat = other.coef_mat;
    curr_zi = other.curr_zi;
    saved_ti = other.saved_ti;
    ai = other.ai;
    coef_n = other.coef_n;
    coef_d = other.coef_d;
    non_solved_degs_den = other.non_solved_degs_den;
    non_solved_degs_num = other.non_solved_degs_num;
    saved_num_num = other.saved_num_num;
    saved_num_den = other.saved_num_den;
    max_deg_num = other.max_deg_num;
    max_deg_den = other.max_deg_den;
    curr_deg_num = other.curr_deg_num;
    curr_deg_den = other.curr_deg_den;
    curr_zi_order_num = other.curr_zi_order_num;
    curr_zi_order_den = other.curr_zi_order_den;
    tmp_solved_coefs_num = other.tmp_solved_coefs_num;
    tmp_solved_coefs_den = other.tmp_solved_coefs_den;
    result = other.result;
    ti = other.ti;
    g_ni = other.g_ni;
    g_di = other.g_di;
    combined_ni = other.combined_ni;
    combined_di = other.combined_di;
    coef_mat_num = other.coef_mat_num;
    coef_mat_den = other.coef_mat_den;
    solved_num = other.solved_num;
    solved_den = other.solved_den;
    solved_degs_num = other.solved_degs_num;
    solved_degs_den = other.solved_degs_den;
    min_deg_den_vec = other.min_deg_den_vec;
    is_singular_system = other.is_singular_system;
    queue = other.queue;
    const_den = other.const_den;
    tag = other.tag;
    sub_num = other.sub_num;
    sub_den = other.sub_den;
    parsed_variables = other.parsed_variables;
    curr_parsed_variable = other.curr_parsed_variable;

    done = other.done;
    new_prime = other.new_prime;
    check = other.check;
    use_chinese_remainder = other.use_chinese_remainder;
    curr_zi_order = other.curr_zi_order;
    prime_number = other.prime_number;
    num_eqn = other.num_eqn;
    n = other.n;
    type = other.type;
    zi = other.zi;
    combined_prime = other.combined_prime;
  }

  RatReconst::RatReconst(RatReconst && other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    first_run = std::move(other.first_run);
    coef_mat = std::move(other.coef_mat);
    curr_zi = std::move(other.curr_zi);
    saved_ti = std::move(other.saved_ti);
    ai = std::move(other.ai);
    coef_n = std::move(other.coef_n);
    coef_d = std::move(other.coef_d);
    non_solved_degs_den = std::move(other.non_solved_degs_den);
    non_solved_degs_num = std::move(other.non_solved_degs_num);
    saved_num_num = std::move(other.saved_num_num);
    saved_num_den = std::move(other.saved_num_den);
    max_deg_num = std::move(other.max_deg_num);
    max_deg_den = std::move(other.max_deg_den);
    curr_deg_num = std::move(other.curr_deg_num);
    curr_deg_den = std::move(other.curr_deg_den);
    curr_zi_order_num = std::move(other.curr_zi_order_num);
    curr_zi_order_den = std::move(other.curr_zi_order_den);
    tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
    tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
    result = std::move(other.result);
    ti = std::move(other.ti);
    g_ni = std::move(other.g_ni);
    g_di = std::move(other.g_di);
    combined_ni = std::move(other.combined_ni);
    combined_di = std::move(other.combined_di);
    coef_mat_num = std::move(other.coef_mat_num);
    coef_mat_den = std::move(other.coef_mat_den);
    solved_num = std::move(other.solved_num);
    solved_den = std::move(other.solved_den);
    solved_degs_num = std::move(other.solved_degs_num);
    solved_degs_den = std::move(other.solved_degs_den);
    min_deg_den_vec = std::move(other.min_deg_den_vec);
    is_singular_system = std::move(other.is_singular_system);
    queue = std::move(other.queue);
    const_den = std::move(other.const_den);
    tag = std::move(other.tag);
    sub_num = std::move(other.sub_num);
    sub_den = std::move(other.sub_den);
    parsed_variables = std::move(other.parsed_variables);
    curr_parsed_variable = std::move(other.curr_parsed_variable);

    done = std::move(other.done);
    new_prime = std::move(other.new_prime);
    check = std::move(other.check);
    use_chinese_remainder = std::move(other.use_chinese_remainder);
    curr_zi_order = std::move(other.curr_zi_order);
    prime_number = std::move(other.prime_number);
    num_eqn = std::move(other.num_eqn);
    n = std::move(other.n);
    type = std::move(other.type);
    zi = std::move(other.zi);
    combined_prime = std::move(other.combined_prime);
  }

  RatReconst& RatReconst::operator=(const RatReconst& other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      first_run = other.first_run;
      coef_mat = other.coef_mat;
      curr_zi = other.curr_zi;
      saved_ti = other.saved_ti;
      ai = other.ai;
      coef_n = other.coef_n;
      coef_d = other.coef_d;
      saved_num_num = other.saved_num_num;
      saved_num_den = other.saved_num_den;
      non_solved_degs_den = other.non_solved_degs_den;
      non_solved_degs_num = other.non_solved_degs_num;
      max_deg_num = other.max_deg_num;
      max_deg_den = other.max_deg_den;
      curr_deg_num = other.curr_deg_num;
      curr_deg_den = other.curr_deg_den;
      curr_zi_order_num = other.curr_zi_order_num;
      curr_zi_order_den = other.curr_zi_order_den;
      tmp_solved_coefs_num = other.tmp_solved_coefs_num;
      tmp_solved_coefs_den = other.tmp_solved_coefs_den;
      result = other.result;
      ti = other.ti;
      g_ni = other.g_ni;
      g_di = other.g_di;
      combined_ni = other.combined_ni;
      combined_di = other.combined_di;
      coef_mat_num = other.coef_mat_num;
      coef_mat_den = other.coef_mat_den;
      solved_num = other.solved_num;
      solved_den = other.solved_den;
      solved_degs_num = other.solved_degs_num;
      solved_degs_den = other.solved_degs_den;
      min_deg_den_vec = other.min_deg_den_vec;
      is_singular_system = other.is_singular_system;
      queue = other.queue;
      const_den = other.const_den;
      tag = other.tag;
      sub_num = other.sub_num;
      sub_den = other.sub_den;
      parsed_variables = other.parsed_variables;
      curr_parsed_variable = other.curr_parsed_variable;

      done = other.done;
      new_prime = other.new_prime;
      check = other.check;
      use_chinese_remainder = other.use_chinese_remainder;
      curr_zi_order = other.curr_zi_order;
      prime_number = other.prime_number;
      num_eqn = other.num_eqn;
      n = other.n;
      type = other.type;
      zi = other.zi;
      combined_prime = other.combined_prime;
    }

    return *this;
  }

  RatReconst& RatReconst::operator=(RatReconst && other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      first_run = std::move(other.first_run);
      coef_mat = std::move(other.coef_mat);
      curr_zi = std::move(other.curr_zi);
      saved_ti = std::move(other.saved_ti);
      ai = std::move(other.ai);
      coef_n = std::move(other.coef_n);
      coef_d = std::move(other.coef_d);
      saved_num_num = std::move(other.saved_num_num);
      saved_num_den = std::move(other.saved_num_den);
      non_solved_degs_den = std::move(other.non_solved_degs_den);
      non_solved_degs_num = std::move(other.non_solved_degs_num);
      max_deg_num = std::move(other.max_deg_num);
      max_deg_den = std::move(other.max_deg_den);
      curr_deg_num = std::move(other.curr_deg_num);
      curr_deg_den = std::move(other.curr_deg_den);
      curr_zi_order_num = std::move(other.curr_zi_order_num);
      curr_zi_order_den = std::move(other.curr_zi_order_den);
      tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
      tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
      result = std::move(other.result);
      ti = std::move(other.ti);
      g_ni = std::move(other.g_ni);
      g_di = std::move(other.g_di);
      combined_ni = std::move(other.combined_ni);
      combined_di = std::move(other.combined_di);
      coef_mat_num = std::move(other.coef_mat_num);
      coef_mat_den = std::move(other.coef_mat_den);
      solved_num = std::move(other.solved_num);
      solved_den = std::move(other.solved_den);
      solved_degs_num = std::move(other.solved_degs_num);
      solved_degs_den = std::move(other.solved_degs_den);
      min_deg_den_vec = std::move(other.min_deg_den_vec);
      is_singular_system = std::move(other.is_singular_system);
      queue = std::move(other.queue);
      const_den = std::move(other.const_den);
      tag = std::move(other.tag);
      sub_num = std::move(other.sub_num);
      sub_den = std::move(other.sub_den);
      parsed_variables = std::move(other.parsed_variables);
      curr_parsed_variable = std::move(other.curr_parsed_variable);

      done = std::move(other.done);
      new_prime = std::move(other.new_prime);
      check = std::move(other.check);
      use_chinese_remainder = std::move(other.use_chinese_remainder);
      curr_zi_order = std::move(other.curr_zi_order);
      prime_number = std::move(other.prime_number);
      num_eqn = std::move(other.num_eqn);
      n = std::move(other.n);
      type = std::move(other.type);
      zi = std::move(other.zi);
      combined_prime = std::move(other.combined_prime);
    }

    return *this;
  }

  void RatReconst::disable_shift() {
    shift = std::vector<FFInt> (n, 0);
  }

  void RatReconst::build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);
    FFInt res = (1 - const_den) * tmp_num;

    for (uint32_t i = 0; i < n; i++) {
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        yis[i] = yis[i] * tmp_ti + shift[i];
      }
    }

    for (int r = 0; r <= max_deg_num; r++) {
      // If the current degree is smaller than the total degree of the polynomial
      // subtract the higher terms to save numerical runs
      if (coef_n.size() > 0 && coef_n[r].is_new_prime())
        res -= coef_n[r].get_result_ff().calc(yis);
      else
        eq.emplace_back(tmp_ti.pow(r));
    }

    for (int r = 1; r <= max_deg_den; r++) {
      if (coef_d.size() > 0 && coef_d[r].is_new_prime()) {
        res += coef_d[r].get_result_ff().calc(yis) * tmp_num;
      } else
        eq.emplace_back(-tmp_ti.pow(r) * tmp_num);
    }

    eq.emplace_back(res);

    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, std::vector<FFInt>& yis) {
    if (!is_singular_system) {
      std::vector<FFInt> eq;
      eq.reserve(num_eqn + 1);

      for (uint32_t i = 0; i < n; i++) {
        yis[i] = yis[i] * tmp_ti + shift[i];
      }

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


      eq.back() += solved_den.calc(yis) * tmp_num;
      eq.back() -= solved_num.calc(yis);

      coef_mat.emplace_back(std::move(eq));
    } else {
      std::vector<FFInt> eq;
      eq.reserve(num_eqn + 1);
      FFInt res = (1 - const_den) * tmp_num;

      for (uint32_t i = 0; i < n; i++) {
        {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          yis[i] = yis[i] * tmp_ti + shift[i];
        }
      }

      for (int r = 0; r <= max_deg_num; r++) {
        // If the current degree is smaller than the total degree of the polynomial
        // subtract the higher terms to save numerical runs
        if (r > curr_deg_num) {
          if (solved_degs_num.find(r) != solved_degs_num.end())
            res -= solved_degs_num[r].calc(yis);
        } else
          eq.emplace_back(tmp_ti.pow(r));
      }

      for (int r = 1; r <= max_deg_den; r++) {
        if (r > curr_deg_den) {
          if (solved_degs_den.find(r) != solved_degs_den.end())
            res += solved_degs_den[r].calc(yis) * tmp_num;
        } else
          eq.emplace_back(-tmp_ti.pow(r) * tmp_num);
      }

      eq.emplace_back(res);

      coef_mat.emplace_back(std::move(eq));
    }
  }

  void RatReconst::generate_anchor_points() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    rand_zi.clear();

    for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 0), 1));
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 1), get_rand()));
    }
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

  void RatReconst::check_for_solved_degs(std::vector<uint32_t>& uni_degs, const bool is_num) {
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

  PolynomialFF RatReconst::solve_transposed_vandermonde(std::vector<std::vector<uint32_t>>& degs,
  const std::vector<std::pair<FFInt, uint32_t>>& nums) {
    uint32_t num_eqn = degs.size();
    std::vector<FFInt> result(num_eqn);

    // calculate base entries of Vandermonde matrix
    std::vector<FFInt> vis;
    vis.reserve(num_eqn);
    std::sort(degs.begin(), degs.end(), std::greater<std::vector<uint32_t>>());

    for (const auto & el : degs) {
      FFInt vi = 1;

      // z_1 is always = 1 which does not matter while determining the coefficient
      for (uint32_t tmp_zi = 2; tmp_zi <= n; tmp_zi++) {
        // curr_zi_ord starts at 1, thus we need to subtract 1 entry
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        vi *= rand_zi[std::make_pair(tmp_zi, 1)].pow(el[tmp_zi - 1]);
      }

      vis.emplace_back(vi);
    }

    // Initialize the coefficient vector of the master polynomial
    std::vector<FFInt> cis(num_eqn);

    // The coefficients of the master polynomial are found by recursion
    // where we have
    // P(Z) = (Z - v_0)*(Z - v_1)*...*(Z - v_{n-1})
    //      =  c_0 + c_1*Z + ... + Z^n
    cis[num_eqn - 1] = -vis[0];

    for (uint32_t i = 1; i < num_eqn; i++) {
      for (uint32_t j = num_eqn - 1 - i; j < num_eqn - 1; j++) {
        cis[j] -= vis[i] * cis[j + 1];
      }

      cis[num_eqn - 1] -= vis[i];
    }

    // Each subfactor in turn is synthetically divided,
    // matrix-multiplied by the right hand-side,
    // and supplied with a denominator (since all vi should be different,
    // there is no additional check if a coefficient in synthetical division
    // leads to a vanishing denominator)
    for (uint32_t i = 0; i < num_eqn; i++) {
      FFInt t = 1;
      FFInt b = 1;
      FFInt s = nums[num_eqn - 1].first;

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += nums[j - 1].first * b;
        t = vis[i] * t + b;
      }

      result[i] = s / t / vis[i];
    }

    // Bring result in canonical form
    ff_map poly;

    for (uint32_t i = 0; i < num_eqn; i ++) {
      poly.emplace(std::make_pair(degs[i], result[i]));
    }

    return PolynomialFF(n, poly);
  }

  FFInt RatReconst::get_rand_zi(uint32_t zi, uint32_t order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return rand_zi.at(std::make_pair(zi, order));
  }

  std::vector<FFInt> RatReconst::get_rand_zi_vec(std::vector<uint32_t> order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    std::vector<FFInt> res {};

    for (uint32_t i = 2; i <= n; i++) {
      res.emplace_back(rand_zi.at(std::make_pair(i, order[i - 2])));
    }

    return res;
  }

  FFInt RatReconst::get_zi_shift(uint32_t zi) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return shift[zi - 1];
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
      for (uint32_t i = 0; i < coef_mat_num[key].size(); i++) {
        auto tmp_pair = coef_mat_num[key][i];

        if (tmp_pair.second < sub_num[curr_deg_num].size()) {
          std::vector<uint32_t> tmp_zi_ord(n - 1, i + 1);
          std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
          yis.emplace(yis.begin(), 1);
          tmp_pair.first -= sub_num[key][tmp_pair.second].calc(yis);
          coef_mat_num[key][i] = tmp_pair;
        }
      }
    }

    solved_degs_num[key] = solve_transposed_vandermonde(non_solved_degs_num[key], coef_mat_num[key]);

    std::vector<uint32_t> zero_deg(n);
    PolynomialFF zero_poly(n, {{zero_deg, 0}});


    for (int i = 0; i < curr_deg_num; i++) {
      if (sub_num[(uint32_t)i].size() == 0)
        sub_num[(uint32_t)i] = {zero_poly};
      else
        sub_num[(uint32_t)i].emplace_back(zero_poly);
    }

    if (curr_deg_num > 0) {
      std::vector<FFInt> tmp_shift;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        tmp_shift = shift;
      }
      PolynomialFF sub_pol = solved_degs_num[key].add_shift(tmp_shift);
      sub_pol -= solved_degs_num[key];

      for (auto & el : sub_pol.coefs) {
        int tmp_deg = 0;

        for (const auto & deg : el.first) tmp_deg += deg;

        if (tmp_deg < curr_deg_num) {
          for (auto & tmp_sub : sub_num[(uint32_t)tmp_deg]) {
            tmp_sub += PolynomialFF(n, {{el.first, el.second}});
          }
        }
      }
    }

    bool found = false;
    uint32_t solved_degs = 1;
    curr_deg_num --;

    if (curr_deg_num > -1) {
      while (!found) {
        if (coef_mat_num.find(curr_deg_num) == coef_mat_num.end()) {
          curr_deg_num --;
          solved_degs ++;
        } else
          found = true;

        if (curr_deg_num == -1)
          found = true;
      }
    }

    tmp_solved_coefs_num += solved_degs;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 1
                - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }
  }

  void RatReconst::set_new_curr_deg_den_singular(uint32_t key) {
    if (curr_deg_den < max_deg_den) {
      for (uint32_t i = 0; i < coef_mat_den[key].size(); i++) {
        auto tmp_pair = coef_mat_den[key][i];

        if (tmp_pair.second < sub_den[curr_deg_den].size()) {
          std::vector<uint32_t> tmp_zi_ord(n - 1, i + 1);
          std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
          yis.emplace(yis.begin(), 1);
          tmp_pair.first -= sub_den[key][tmp_pair.second].calc(yis);
          coef_mat_den[key][i] = tmp_pair;
        }
      }
    }

    solved_degs_den[key] = solve_transposed_vandermonde(non_solved_degs_den[key], coef_mat_den[key]);

    std::vector<uint32_t> zero_deg(n);
    PolynomialFF zero_poly(n, {{zero_deg, 0}});


    for (int i = 0; i < curr_deg_den; i++) {
      if (sub_den[(uint32_t)i].size() == 0)
        sub_den[(uint32_t)i] = {zero_poly};
      else
        sub_den[(uint32_t)i].emplace_back(zero_poly);
    }

    if (curr_deg_den > 0) {
      std::vector<FFInt> tmp_shift;
      {
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        tmp_shift = shift;
      }
      PolynomialFF sub_pol = solved_degs_den[key].add_shift(tmp_shift);
      sub_pol -= solved_degs_den[key];

      for (auto & el : sub_pol.coefs) {
        int tmp_deg = 0;

        for (const auto & deg : el.first) tmp_deg += deg;

        if (tmp_deg < curr_deg_den) {
          for (auto & tmp_sub : sub_den[(uint32_t)tmp_deg]) {
            tmp_sub += PolynomialFF(n, {{el.first, el.second}});
          }
        }
      }
    }

    std::vector<FFInt> tmp_yis(n, 0);
    const_den += sub_den[0].back().calc(tmp_yis);

    bool found = false;
    uint32_t solved_degs = 1;
    curr_deg_den --;

    if (curr_deg_den > 0) {
      while (!found) {
        if (coef_mat_den.find(curr_deg_den) == coef_mat_den.end()) {
          curr_deg_den --;
          solved_degs ++;

        } else
          found = true;

        if (curr_deg_den == 0) {
          curr_deg_den = -1;
          found = true;
        }
      }
    }

    if (curr_deg_den == 0) {
      curr_deg_den = -1;
      found = true;
    }

    tmp_solved_coefs_den += solved_degs;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 1
                - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }
  }

  void RatReconst::set_tag(std::string tag_) {
    tag = tag_;
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
          if (line == "max_deg_num") {
            curr_parsed_variable = MAX_DEG_NUM;
            parsed_variables[MAX_DEG_NUM] = true;
          } else if (line == "max_deg_den") {
            curr_parsed_variable = MAX_DEG_DEN;
            parsed_variables[MAX_DEG_DEN] = true;
          } else if (line == "need_prime_shift") {
            curr_parsed_variable = NEED_PRIME_SHIFT;
            parsed_variables[NEED_PRIME_SHIFT] = true;
          } else if (line == "min_deg_den_vec") {
            curr_parsed_variable = MIN_DEG_DEN_VEC;
            parsed_variables[MIN_DEG_DEN_VEC] = true;
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
          } else {
            switch (curr_parsed_variable) {
              case COMBINED_PRIME: {
                combined_prime = mpz_class(line);
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
                need_prime_shift = std::stoi(line);
                break;
              }

              case MIN_DEG_DEN_VEC: {
                min_deg_den_vec = parse_vector(line);
                n = min_deg_den_vec.size();
                break;
              }

              case G_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'min_deg_den_vec' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                std::vector<mpz_class> tmp_rn = parse_rational_number(line);

                g_ni.emplace(std::make_pair(tmp_vec, RationalNumber(tmp_rn[0], tmp_rn[1])));

                break;
              }

              case G_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'min_deg_den_vec' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                std::vector<mpz_class> tmp_rn = parse_rational_number(line);

                g_di.emplace(std::make_pair(tmp_vec, RationalNumber(tmp_rn[0], tmp_rn[1])));
                break;

              }

              case COMBINED_NI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'min_deg_den_vec' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                combined_ni.emplace(std::make_pair(tmp_vec, mpz_class(line)));

                break;
              }

              case COMBINED_DI: {
                if (n == 0) {
                  ERROR_MSG("Input file is in the wrong order! Need to parse 'min_deg_den_vec' first.");
                  std::exit(-1);
                }

                std::vector<uint32_t> tmp_vec = parse_vector(line, n);
                combined_di.emplace(std::make_pair(tmp_vec, mpz_class(line)));

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
      is_singular_system = need_prime_shift;

      new_prime = true;
      std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
      new_prime = true;

      if (is_singular_system) {
        tmp_solved_coefs_den = 0;
        tmp_solved_coefs_num = 0;
        {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          need_prime_shift = true;
        }

        for (const auto & el : g_ni) add_non_solved_num(el.first);

        for (const auto & el : g_di) add_non_solved_den(el.first);

        curr_deg_num = max_deg_num;
        curr_deg_den = max_deg_den;

        std::unique_lock<std::mutex> lock(mutex_status);
        num_eqn = max_deg_den + max_deg_num + 1;
      } else {
        std::unique_lock<std::mutex> lock(mutex_status);
        num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
      }

      for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

      for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};
    } else {
      ERROR_MSG(std::string("The file '") + file_name + std::string("' could not be found!"));
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

  std::vector< mpz_class > RatReconst::parse_rational_number(std::string& line) {
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

    for (const auto & el : g_ni) add_non_solved_num(el.first);

    for (const auto & el : g_di) add_non_solved_den(el.first);

    curr_deg_num = max_deg_num;
    curr_deg_den = max_deg_den;

    std::unique_lock<std::mutex> lock(mutex_status);
    num_eqn = max_deg_den + max_deg_num + 1;

    for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};

    for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<std::pair<FFInt, uint32_t>> {};
  }
}
