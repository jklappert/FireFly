#include "Logger.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::need_prime_shift = false;
  bool RatReconst::set_singular_system = false;
  ff_pair_map RatReconst::rand_zi;
  std::mutex RatReconst::mutex_statics;

  RatReconst::RatReconst(uint n_) {
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
          for (auto & el : shift) el = FFInt(std::rand() % 1000000) + FFInt(1);

          curr_zi_order_num = std::vector<uint> (n - 1, 1);
          curr_zi_order_den = std::vector<uint> (n - 1, 1);
        }
      }
    }

    if (n > 1) {
      curr_zi_order = std::vector<uint> (n - 1, 1);
      lock_status.unlock();
      // add a zero to both polynomials to do arithmetic
      ff_map zero_deg {};
      zero_deg.emplace(std::make_pair(std::vector<uint> (n), 0));
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

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime) {
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

  void RatReconst::interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord) {
    if (!done) {

      // Compare if the food is the expected food; if not, store it for later use
      if (feed_zi_ord == curr_zi_order) {
        {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);

          if (!is_singular_system && set_singular_system) is_singular_system = true;
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
              coef_n.clear();
              coef_d.clear();
              combined_di.clear();
              combined_ni.clear();
              combined_prime = 0;
              num_eqn = 0;
              new_prime = false;
              solved_den.coefs.clear();
              solved_num.coefs.clear();
              curr_zi_order.clear();
              saved_num_num.clear();
              saved_num_den.clear();
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
          const uint i = ti.size() - 1;

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
              for (uint i = 0; i < curr_zi_order.size(); i ++) {
                std::unique_lock<std::mutex> lock_statics(mutex_statics);
                yis.emplace_back(rand_zi[std::make_pair(i + 2, curr_zi_order[i])]);

                if (prime_number > 0 && !is_singular_system)
                  yis.back() *= tmp_ti;
              }
            }

            yis.insert(yis.begin(), FFInt(1));

            if (prime_number > 0  && !is_singular_system)
              yis.front() *= tmp_ti;

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

              for (uint j = 0; j < n; j++) {
                shift[j] = FFInt(std::rand() % 1000000) + FFInt(1);
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
            std::vector<uint> zero_deg(n);
            ff_map zero_mon = {{zero_deg, 0}};

            // save the current results to the map to access them later
            for (int i = 0; i <= (int)(max_deg_num - tmp_solved_coefs_num); i++) {
              if (first_run) {
                {
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  PolyReconst rec(n - 1, (uint) i, true);

                  if (rec.is_rand_zi_empty()) {
                    std::vector<FFInt> anchor_points {};

                    for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
                      anchor_points.emplace_back(rand_zi[std::make_pair(tmp_zi, 1)]);
                    }

                    rec.set_anchor_points(anchor_points);
                  }

                  coef_n.emplace(std::make_pair((uint) i, std::move(rec)));
                }

                sub_num.emplace(std::make_pair(i, PolynomialFF(n, zero_mon)));
              }

              if (i <= curr_deg_num) {
                // this saves some memory since we only need one numerical value
                // for the constant coefficient
                if (i == 0 && first_run) {
                  std::vector<uint> key = {(uint) i, zi};
                  saved_num_num[curr_zi_order][key] = num_coef[ {(uint) i}];
                } else {
                  std::vector<uint> key = {(uint) i, zi};
                  saved_num_num[curr_zi_order][key] = num_coef[ {(uint) i}];
                }
              }
            }

            for (uint i = 1; i <= max_deg_den - tmp_solved_coefs_den; i++) {
              if (first_run) {
                {
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  PolyReconst rec(n - 1, i, true);

                  if (rec.is_rand_zi_empty()) {
                    std::vector<FFInt> anchor_points {};

                    for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi ++)
                      anchor_points.emplace_back(rand_zi[std::make_pair(tmp_zi, 1)]);

                    rec.set_anchor_points(anchor_points);
                  }

                  coef_d.emplace(std::make_pair(i, std::move(rec)));
                }

                sub_den.emplace(std::make_pair(i, PolynomialFF(n, zero_mon)));

                if (i == 1)
                  sub_den.emplace(std::make_pair(0, PolynomialFF(n, zero_mon)));
              }

              if ((int) i <= curr_deg_den) {
                std::vector<uint> key = {i, zi};
                saved_num_den[curr_zi_order][key] = den_coef[ {i}];
              }
            }

            if (first_run) first_run = false;

            uint zi_num = 0;

            if (curr_deg_num > 0) zi_num = coef_n[curr_deg_num].get_zi() + 1;

            uint zi_den = 0;

            if (curr_deg_den > 0) zi_den = coef_d[curr_deg_den].get_zi() + 1;

            // reconstruct the numerator
            if (curr_deg_num >= 0) {
              PolyReconst rec = coef_n[curr_deg_num];

              if (rec.get_zi_order() == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end())) {
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

              if (rec.get_zi_order() == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end())) {
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

                  if (rec_new.get_zi_order() == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end())) {
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
              std::vector<uint> zi_order_num_rev = curr_zi_order_num;
              std::reverse(zi_order_num_rev.begin(), zi_order_num_rev.end());
              std::vector<uint> zi_order_den_rev = curr_zi_order_den;
              std::reverse(zi_order_den_rev.begin(), zi_order_den_rev.end());

              if (zi_order_num_rev > zi_order_den_rev) {
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
                PolynomialFF res = el.second.get_result_ff().homogenize(el.first);

                if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  numerator += res;
              }

              for (auto & el : coef_d) {
                PolynomialFF res = el.second.get_result_ff().homogenize(el.first);

                if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0))
                  denominator += res;
              }

              FFInt terminator = 0;
              // check if one can normalize to a single univariate degree, if so
              // set the corresponding degree in equalizer_degree. Check constants
              // first and proceed with denominator/numerator.
              // if the denominator is just a constant, there is no corresponding
              // PolyReconst object. Thus we set the minimal degree to a zero tuple
              FFInt const_shift = sub_den[0].calc(std::vector<FFInt> (n, 0));

              if (const_shift != 1) {
                ff_map dummy_map {};
                terminator = FFInt(1) - const_shift;
                dummy_map.emplace(std::make_pair(std::vector<uint> (n, 0), terminator));
                denominator = denominator + PolynomialFF(n, dummy_map);
              } else if (numerator.coefs.find(std::vector<uint> (n, 0)) != numerator.coefs.end()) {
                terminator = numerator.coefs[std::vector<uint> (n, 0)];
              } else {
                for (const auto & el : denominator.coefs) {
                  std::vector<uint> degs_reverse = el.first;
                  std::reverse(degs_reverse.begin(), degs_reverse.end());

                  add_non_solved_den(el.first);

                  if (min_deg_den_vec.empty())
                    min_deg_den_vec = el.first;
                  else {
                    std::reverse(min_deg_den_vec.begin(), min_deg_den_vec.end());

                    if (min_deg_den_vec > degs_reverse) {
                      min_deg_den_vec = degs_reverse;
                    }

                    std::reverse(min_deg_den_vec.begin(), min_deg_den_vec.end());
                  }
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
              for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
                auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 2]);
                std::unique_lock<std::mutex> lock_statics(mutex_statics);
                rand_zi.emplace(std::make_pair(key, rand_zi[std::make_pair(tmp_zi, 1)].pow(key.second)));
              }
            }

            return;
          } else {
            // build multivariate Vandermonde systems and evaluate them if possible
            for (const auto & sol : canonical.first) {
              uint key = sol.first[0];

              if (coef_mat_num.find(key) != coef_mat_num.end()) {
                if (coef_mat_num[key].size() < non_solved_degs_num[key].size())
                  coef_mat_num[key].emplace_back(sol.second);

                // Solve multivariate gauss system for corresponding degree,
                // remove entry from non_solved_degs and add it to solve_degs
                if ((key == curr_deg_num || !is_singular_system) && coef_mat_num[key].size() == non_solved_degs_num[key].size()) {
                  if (!is_singular_system) {
                    solved_num += solve_transposed_vandermonde(non_solved_degs_num[key], coef_mat_num[key]);
                    non_solved_degs_num.erase(key);
                  } else {
                    set_new_curr_deg_num_singular(key);
                    bool cannot_solve = false;

                    while (!cannot_solve) {
                      std::vector<uint> deleted_degs {};

                      for (const auto & mat : coef_mat_num) {
                        uint tmp_key = mat.first;

                        if (tmp_key == curr_deg_num && mat.second.size() == non_solved_degs_num[tmp_key].size()) {
                          set_new_curr_deg_num_singular(tmp_key);
                          deleted_degs.emplace_back(tmp_key);
                        }
                      }

                      for (const auto & el : deleted_degs) {
                        coef_mat_num.erase(el);
                      }

                      if (deleted_degs.empty()) cannot_solve = true;
                    }
                  }

                  coef_mat_num.erase(key);
                }
              }
            }

            for (const auto & sol : canonical.second) {
              uint key = sol.first[0];

              if (coef_mat_den.find(key) != coef_mat_den.end()) {
                if (coef_mat_den[key].size() < non_solved_degs_den[key].size())
                  coef_mat_den[key].emplace_back(sol.second);

                // Solve multivariate Vandermonde system for corresponding degree,
                // remove entry from non_solved_degs and add it to solve_degs
                if ((key == curr_deg_den || !is_singular_system) && coef_mat_den[key].size() == non_solved_degs_den[key].size()) {
                  if (!is_singular_system) {
                    solved_den += solve_transposed_vandermonde(non_solved_degs_den[key], coef_mat_den[key]);
                    non_solved_degs_den.erase(key);
                  } else {
                    set_new_curr_deg_den_singular(key);
                    bool cannot_solve = false;

                    while (!cannot_solve) {
                      std::vector<uint> deleted_degs {};

                      for (const auto & mat : coef_mat_den) {
                        uint tmp_key = mat.first;

                        if (tmp_key == curr_deg_den && mat.second.size() == non_solved_degs_den[tmp_key].size()) {
                          set_new_curr_deg_den_singular(tmp_key);
                          deleted_degs.emplace_back(tmp_key);
                        }
                      }

                      for (const auto & el : deleted_degs) {
                        coef_mat_den.erase(el);
                      }

                      if (deleted_degs.empty()) cannot_solve = true;
                    }
                  }

                  coef_mat_den.erase(key);
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

              solved_num.coefs.erase(std::vector<uint> (n, 0));
              solved_den.coefs.erase(std::vector<uint> (n, 0));
//              std::cout << solved_num << solved_den;

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
              zero_deg.emplace(std::make_pair(std::vector<uint> (n), 0));
              solved_degs_num.clear();
              solved_degs_den.clear();
              solved_num.coefs = zero_deg;
              solved_den.coefs = zero_deg;
            } else {
              // increase zi order by 1
              {
                std::unique_lock<std::mutex> lock(mutex_status);
                std::transform(curr_zi_order.begin(), curr_zi_order.end(),
                curr_zi_order.begin(), [](uint x) {return x + 1;});

                for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
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

  std::tuple<int, uint, std::vector<uint>> RatReconst::feed_poly(int curr_deg,
                                                                 uint max_deg, std::unordered_map<uint, PolyReconst>& coef,
                                                                 PolyReconst& rec, ff_map_map& saved_num,
  std::unordered_map<uint, PolynomialFF>& sub_save, bool is_num) {
    uint tmp_zi = rec.get_zi() + 1;
    std::vector<uint> tmp_zi_ord = curr_zi_order;

    while (!rec.is_new_prime()) {
      try {
        std::vector<uint> key = {(uint) curr_deg, tmp_zi};
        FFInt food = saved_num.at(tmp_zi_ord).at(key);
        // delete unused saved data
        saved_num[tmp_zi_ord].erase(key);
        // set random values for the yis
        std::vector<FFInt> yis {};

        for (uint i = 0; i < tmp_zi_ord.size(); i ++) {
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          yis.emplace_back(rand_zi[std::make_pair(i + 2, tmp_zi_ord[i])]);
        }

        // feed to PolyReconst
        // since the constant is just a constant, we do not have to get mutliple
        // numerical values to reconstruct the coefficient
        if (curr_deg == 0) {
          while (!rec.is_new_prime()) {
            if (curr_deg == (int) max_deg)
              rec.feed(yis, food);
            else {
              yis.emplace(yis.begin(), FFInt(1));
              FFInt num_subtraction = sub_save[curr_deg].calc(yis);
              yis.erase(yis.begin());
              rec.feed(yis, food - num_subtraction);
            }
          }
        } else {
          if (curr_deg == (int) max_deg)
            rec.feed(yis, food);
          else {
            yis.emplace(yis.begin(), FFInt(1));
            FFInt num_subtraction = sub_save[curr_deg].calc(yis);
            yis.erase(yis.begin());

            rec.feed(yis, food - num_subtraction);
          }

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

          // check if the polynomial is zero which we then can omit for further
          // calculations
          if (!(res.coefs.size() == 1 && res.coefs.begin()->second == 0)) {
            PolynomialFF sub_pol = rec.get_result_ff().homogenize(curr_deg).add_shift(shift);
            sub_pol -= rec.get_result_ff().homogenize(curr_deg);

            for (auto & el : sub_pol.coefs) {
              uint tmp_deg = 0;

              for (auto & n : el.first) {
                tmp_deg += n;
              }

              ff_map monomial = {{el.first, el.second}};
              sub_save[tmp_deg] += PolynomialFF(n, monomial);
            }
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
    std::vector<uint> tmp_deg_num {};
    std::vector<uint> tmp_deg_den {};

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

    if (!use_chinese_remainder) {
      combined_ni = tmp.first;
      combined_di = tmp.second;

      // if the coefficient is not a rational number thus divided by 1,
      // it will not change in the next run and can be omitted to save
      // numerical runs
      mpz_map combined_ni_back = combined_ni;

      for (auto & c_ni : combined_ni_back) {
        try {
          RationalNumber rn = get_rational_coef(c_ni.second, combined_prime);

          if (rn.numerator == c_ni.second && rn.denominator == 1)
            remove_ni(c_ni.first, rn);
          else
            add_non_solved_num(c_ni.first);
        } catch (std::exception& e) {
          add_non_solved_num(c_ni.first);
        }
      }

      mpz_map combined_di_back = combined_di;

      for (auto & c_di : combined_di_back) {
        try {
          RationalNumber rn = get_rational_coef(c_di.second, combined_prime);

          if (rn.numerator == c_di.second && rn.denominator == 1)
            remove_di(c_di.first, rn);
          else
            add_non_solved_den(c_di.first);
        } catch (std::exception& e) {
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
          try {
            RationalNumber last_rn = get_rational_coef(c_ni.second, combined_prime_back);
            RationalNumber curr_rn = get_rational_coef(combined_ni[c_ni.first], combined_prime);

            if (last_rn == curr_rn)
              remove_ni(c_ni.first, curr_rn);
            else
              add_non_solved_num(c_ni.first);
          } catch (std::exception& e) {

            if (c_ni.second == combined_ni[c_ni.first]) {
              RationalNumber rn = RationalNumber(c_ni.second, 1);
              remove_ni(c_ni.first, rn);
            } else
              add_non_solved_num(c_ni.first);
          }
        }
      }

      for (auto & c_di : combined_di_back) {
        if (g_di.find(c_di.first) == g_di.end()) {
          try {
            RationalNumber last_rn = get_rational_coef(c_di.second, combined_prime_back);
            RationalNumber curr_rn = get_rational_coef(combined_di[c_di.first], combined_prime);

            if (last_rn == curr_rn)
              remove_di(c_di.first, curr_rn);
            else
              add_non_solved_den(c_di.first);
          } catch (std::exception& e) {

            if (c_di.second == combined_di[c_di.first]) {
              RationalNumber rn = RationalNumber(c_di.second, 1);

              remove_di(c_di.first, rn);
            } else
              add_non_solved_den(c_di.first);
          }
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

    sub_num.clear();
    sub_den.clear();

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

      std::vector<uint> zero_deg(n);
      ff_map zero_mon = {{zero_deg, 0}};
      PolynomialFF zero_poly(n, zero_mon);

      for (int i = 0; i <= std::max(max_deg_num, max_deg_den); i++) {
        if (i <= max_deg_num)
          sub_num.emplace(std::make_pair(i, zero_poly));

        if (i >= 1 && i <= max_deg_den) {
          sub_den.emplace(std::make_pair(i, zero_poly));

          if (i == 1)
            sub_den.emplace(std::make_pair(0, zero_poly));
        }
      }

      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 1
                - tmp_solved_coefs_num - tmp_solved_coefs_den;
    } else {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
    }

    for (const auto & el : non_solved_degs_num) coef_mat_num[el.first] = std::vector<FFInt> {};

    for (const auto & el : non_solved_degs_den) coef_mat_den[el.first] = std::vector<FFInt> {};
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
      throw std::runtime_error("Trying to access unfinished result.");
    }
  }

  bool RatReconst::rec_rat_coef() {
    bool run_test = true;
    std::vector<const std::vector<uint>*> promoted_n;
    std::vector<const std::vector<uint>*> promoted_d;

    for (const auto & ci : combined_ni) {
      mpz_class a = ci.second;

      try {
        g_ni[ci.first] = get_rational_coef(a, combined_prime);
        promoted_n.emplace_back(&ci.first);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    for (const auto & ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di[ci.first] = get_rational_coef(a, combined_prime);
        promoted_d.emplace_back(&ci.first);
      } catch (const std::exception&) {
        run_test = false;
        break;
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

  FFInt RatReconst::comp_fyi(uint i, uint ip, const FFInt& y) {
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

    const std::vector<uint> min_power = {0};
    denominator.emplace(std::make_pair(std::move(min_power), FFInt(1)));

    int terms_num = max_deg_num - tmp_solved_coefs_num;

    if (terms_num == -1) {
      numerator.emplace(std::make_pair(std::move(min_power), FFInt(1)));
    } else {
      for (int i = 0; i <= terms_num; i ++) {
        std::vector<uint> power = { (uint) i};
        numerator.emplace(std::make_pair(std::move(power), results[i]));
      }
    }

    for (uint i = 1; i <= max_deg_den - tmp_solved_coefs_den; i ++) {
      std::vector<uint> power = {i};
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

    uint counter = 0;

    if (!is_singular_system) {
      for (const auto & el : non_solved_degs_num) {
        std::vector<uint> power = {el.first};
        numerator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }

      for (const auto & el : non_solved_degs_den) {
        std::vector<uint> power = {el.first};
        denominator.emplace(std::make_pair(std::move(power), results[counter]));
        counter ++;
      }
    } else {
      int terms_num = max_deg_num - tmp_solved_coefs_num;

      if (terms_num > 0) {
        for (uint i = 0; i <= terms_num; i++) {
          if (non_solved_degs_num.find(i) != non_solved_degs_num.end()) {
            std::vector<uint> power = {i};
            numerator.emplace(std::make_pair(std::move(power), results[counter]));
          }

          counter ++;
        }
      }

      for (uint i = 1; i <= max_deg_den - tmp_solved_coefs_den; i++) {
        if (non_solved_degs_den.find(i) != non_solved_degs_den.end()) {
          std::vector<uint> power = {i};
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
      std::vector<uint> zero_deg = {0};
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

  bool RatReconst::test_guess(const FFInt& num) {
    ff_map g_ff_ni = convert_to_ffint(g_ni);
    ff_map g_ff_di = convert_to_ffint(g_di);
    PolynomialFF g_ny(n, g_ff_ni);
    PolynomialFF g_dy(n, g_ff_di);
    std::vector<FFInt> yis = std::vector<FFInt> (n);
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yis[0] = ti[0] + shift[0];

      for (uint i = 1; i < n; i++) {
        yis[i] = ti[0] * rand_zi[std::make_pair(i + 1, curr_zi_order[i - 1])] + shift[i];
      }
    }

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;
  }

  void RatReconst::remove_ni(const std::vector<uint>& deg_vec, RationalNumber& rn) {
    g_ni[deg_vec] =  rn;
    combined_ni.erase(deg_vec);
  }

  void RatReconst::remove_di(const std::vector<uint>& deg_vec, RationalNumber& rn) {
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
    sub_num = other.sub_num;
    sub_den = other.sub_den;
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
    sub_num = std::move(other.sub_num);
    sub_den = std::move(other.sub_den);
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
      sub_num = other.sub_num;
      sub_den = other.sub_den;
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
      sub_num = std::move(other.sub_num);
      sub_den = std::move(other.sub_den);
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

  void RatReconst::build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);
    FFInt res = tmp_num;
    std::vector<FFInt> yis_wo_t = yis;
    yis_wo_t.erase(yis_wo_t.begin());

    for (int r = 0; r <= max_deg_num; r++) {
      // If the current degree is smaller than the total degree of the polynomial
      // subtract the higher terms to save numerical runs
      if (shift[0] != 0 && coef_n.size() > 0 && coef_n[r].is_new_prime())
        res -= (coef_n[r].get_result_ff().calc(yis_wo_t) + sub_num[r].calc(yis)) * tmp_ti.pow(r);
      else
        eq.emplace_back(tmp_ti.pow(r));
    }

    for (int rp = 1; rp <= max_deg_den; rp++) {
      // If the current degree is smaller than the total degree of the polynomial
      // subtract the higher terms to save numerical runs
      if (shift[0] != 0  && coef_d.size() > 0 && coef_d[rp].is_new_prime())
        res += (coef_d[rp].get_result_ff().calc(yis_wo_t) + sub_den[rp].calc(yis)) * tmp_ti.pow(rp) * tmp_num;
      else
        eq.emplace_back(-tmp_ti.pow(rp) * tmp_num);
    }

    // The lowest degree in univariate Gauss for multivariate polynomial
    // reconstruction is always zero. Hence, we just need to emplace the
    // evaluation of f(yis) at this point
    eq.emplace_back(res);

    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    if (!is_singular_system) {
      std::vector<FFInt> eq;
      eq.reserve(num_eqn + 1);

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

      for (const auto & el : g_ni) {
        mpz_class tmp_1 = el.second.numerator % FFInt::p;
        mpz_class tmp_2 = el.second.denominator % FFInt::p;

        if (tmp_1 < 0) tmp_1 = FFInt::p + tmp_1;

        FFInt coef = FFInt(std::stoull(tmp_1.get_str())) / FFInt(std::stoull(tmp_2.get_str()));

        for (uint i = 0; i < n; i++) {
          coef *= yis[i].pow(el.first[i]);
        }

        eq.back() -= coef;
      }

      FFInt sol_den = 0;

      for (const auto & el : g_di) {
        mpz_class tmp_1 = el.second.numerator % FFInt::p;
        mpz_class tmp_2 = el.second.denominator % FFInt::p;

        if (tmp_1 < 0) tmp_1 = FFInt::p + tmp_1;

        FFInt tmp_coef = FFInt(std::stoull(tmp_1.get_str())) / FFInt(std::stoull(tmp_2.get_str()));

        for (uint i = 0; i < n; i++) {
          tmp_coef *= yis[i].pow(el.first[i]);
        }

        sol_den += tmp_coef;
      }

      if (n > 1) {
        sol_den += solved_den.calc(yis);
        eq.back() -= solved_num.calc(yis);
      }

      sol_den *= tmp_num;

      eq.back() += sol_den;

      coef_mat.emplace_back(std::move(eq));
    } else {
      std::vector<FFInt> eq;
      eq.reserve(num_eqn + 1);
      FFInt res = tmp_num;

      for (int r = 0; r <= max_deg_num; r++) {
        // If the current degree is smaller than the total degree of the polynomial
        // subtract the higher terms to save numerical runs
        if (r > curr_deg_num)
          res -= (solved_degs_num[r].calc(yis) + sub_num[r].calc(yis)) * tmp_ti.pow(r);
        else
          eq.emplace_back(tmp_ti.pow(r));
      }

      for (int r = 1; r <= max_deg_den; r++) {
        // If the current degree is smaller than the total degree of the polynomial
        // subtract the higher terms to save numerical runs
        if (r > curr_deg_den)
          res += tmp_num * (solved_degs_den[r].calc(yis) + sub_den[r].calc(yis)) * tmp_ti.pow(r);
        else
          eq.emplace_back(-tmp_ti.pow(r) * tmp_num);
      }

      eq.emplace_back(res);

      coef_mat.emplace_back(std::move(eq));
    }
  }

  void RatReconst::generate_anchor_points() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    rand_zi.clear();

    for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi ++) {
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 0), 1));
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 1), get_rand()));
    }
  }

  void RatReconst::add_non_solved_num(const std::vector<uint>& deg) {
    uint degree = 0;

    for (const auto & el : deg) degree += el;

    non_solved_degs_num[degree].emplace_back(deg);
  }

  void RatReconst::add_non_solved_den(const std::vector<uint>& deg) {
    uint degree = 0;

    for (const auto & el : deg) degree += el;

    non_solved_degs_den[degree].emplace_back(deg);
  }

  void RatReconst::check_for_solved_degs(std::vector<uint>& uni_degs, const bool is_num) {
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

  PolynomialFF RatReconst::solve_transposed_vandermonde(std::vector<std::vector<uint>>& degs,
  const std::vector<FFInt>& nums) {
    uint num_eqn = degs.size();
    std::vector<FFInt> result(num_eqn);

    if (num_eqn == 1) {
      FFInt vi = 1;

      for (const auto & el : degs) {
        for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi++) {
          // curr_zi_ord starts at 1, thus we need to subtract 1 entry
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          vi *= rand_zi[std::make_pair(tmp_zi, 1)].pow(el[tmp_zi - 1]);
        }
      }

      result[0] = nums[0] * 1 / vi;
    } else {
      // calculate base entries of Vandermonde matrix
      std::vector<FFInt> vis;
      vis.reserve(num_eqn);
      std::sort(degs.begin(), degs.end(), std::greater<std::vector<uint>>());

      for (const auto & el : degs) {
        FFInt vi = 1;

        // z_1 is always = 1 which does not matter while determining the coefficient
        for (uint tmp_zi = 2; tmp_zi <= n; tmp_zi++) {
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

      for (uint i = 1; i < num_eqn; i++) {
        for (uint j = num_eqn - 1 - i; j < num_eqn - 1; j++) {
          cis[j] -= vis[i] * cis[j + 1];
        }

        cis[num_eqn - 1] -= vis[i];
      }

      // Each subfactor in turn is synthetically divided,
      // matrix-multiplied by the right hand-side,
      // and supplied with a denominator (since all vi should be different,
      // there is no additional check if a coefficient in synthetical division
      // leads to a vanishing denominator)
      for (uint i = 0; i < num_eqn; i++) {
        FFInt t = 1;
        FFInt b = 1;
        FFInt s = nums[num_eqn - 1];

        for (int j = num_eqn - 1; j > 0; j--) {
          b = cis[j] + vis[i] * b;
          s += nums[j - 1] * b;
          t = vis[i] * t + b;
        }

        result[i] = s / t / vis[i];
      }
    }

    // Bring result in canonical form
    ff_map poly;

    for (uint i = 0; i < num_eqn; i ++) {
      poly.emplace(std::make_pair(degs[i], result[i]));
    }

    return PolynomialFF(n, poly);
  }

  FFInt RatReconst::get_rand_zi(uint zi, uint order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return rand_zi.at(std::make_pair(zi, order));
  }

  std::vector<FFInt> RatReconst::get_rand_zi_vec(std::vector<uint> order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    std::vector<FFInt> res {};

    for (uint i = 2; i <= n; i++) {
      res.emplace_back(rand_zi.at(std::make_pair(i, order[i - 2])));
    }

    return res;
  }

  FFInt RatReconst::get_zi_shift(uint zi) {
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

  void RatReconst::set_new_curr_deg_num_singular(uint key) {
    if (curr_deg_num < max_deg_num) {
      for (uint i = 0; i < coef_mat_num[key].size(); i++) {
        std::vector<uint> tmp_zi_ord(n - 1, i + 1);
        std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
        yis.insert(yis.begin(), 1);
        coef_mat_num[key][i] -= sub_num[key].calc(yis);
      }
    }

    solved_degs_num[key] = solve_transposed_vandermonde(non_solved_degs_num[key], coef_mat_num[key]);

    if (curr_deg_num > 0) {
      PolynomialFF sub_pol = solved_degs_num[key].add_shift(shift);
      sub_pol -= solved_degs_num[key];

      for (auto & el : sub_pol.coefs) {
        uint tmp_deg = 0;

        for (auto & n : el.first) {
          tmp_deg += n;
        }

        ff_map monomial = {{el.first, el.second}};
        sub_num[tmp_deg] += PolynomialFF(n, monomial);
      }
    }

    bool found = false;
    uint solved_degs = 1;
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

  void RatReconst::set_new_curr_deg_den_singular(uint key) {
    if (curr_deg_den < max_deg_den) {
      for (uint i = 0; i < coef_mat_den[key].size(); i++) {
        std::vector<uint> tmp_zi_ord(n - 1, i + 1);
        std::vector<FFInt> yis = get_rand_zi_vec(tmp_zi_ord);
        yis.insert(yis.begin(), 1);
        coef_mat_den[key][i] -= sub_den[key].calc(yis);
      }
    }

    solved_degs_den[key] = solve_transposed_vandermonde(non_solved_degs_den[key], coef_mat_den[key]);

    if (curr_deg_den > 0) {
      PolynomialFF sub_pol = solved_degs_den[key].add_shift(shift);
      sub_pol -= solved_degs_den[key];

      for (auto & el : sub_pol.coefs) {
        uint tmp_deg = 0;

        for (auto & n : el.first) {
          tmp_deg += n;
        }

        ff_map monomial = {{el.first, el.second}};
        sub_den[tmp_deg] += PolynomialFF(n, monomial);
      }
    }

    bool found = false;
    uint solved_degs = 1;
    curr_deg_den --;

    if (curr_deg_den > 0) {
      while (!found) {
        if (coef_mat_den.find(curr_deg_den) == coef_mat_den.end()) {
          curr_deg_den --;
          solved_degs ++;

          if (curr_deg_den == 0) {
            curr_deg_den = -1;
            found = true;
          }
        } else
          found = true;
      }
    }

    tmp_solved_coefs_den += solved_degs;
    {
      std::unique_lock<std::mutex> lock(mutex_status);
      num_eqn = max_deg_den + max_deg_num + 1
                - tmp_solved_coefs_num - tmp_solved_coefs_den;
    }
  }
}



