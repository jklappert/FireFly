#include "Logger.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::shifted = false;
  ff_pair_map RatReconst::rand_zi;
  std::mutex RatReconst::mutex_statics;
  std::vector<FFInt> RatReconst::anchor_points {};

  RatReconst::RatReconst(uint n_) {
    n = n_;
    type = RAT;
    std::unique_lock<std::mutex> lock_status(mutex_status);

    ti.reserve(300);
    ai.reserve(300);
    combined_prime = FFInt::p;

    //std::srand(std::time(nullptr));
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      if (!shifted) {
        shift = std::vector<FFInt> (n);

        if (n > 1) {
          for (auto & el : shift) el = FFInt(std::rand() % 1000000) + FFInt(1);

          curr_zi_order_num = std::vector<uint> (n - 1, 1);
          curr_zi_order_den = std::vector<uint> (n - 1, 1);
          shifted = true;
        }
      }
    }

    if (n > 1) {
      curr_zi_order = std::vector<uint> (n - 1, 1);

      // fill in the rand_vars for zi_order = 1
      if (rand_zi.empty()) {
        generate_anchor_points();
      }
    }
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (fed_prime == prime_number)
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
        lock.lock();
      }
    }

    is_interpolating = false;
  }

  void RatReconst::interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord) {
    // change later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (prime_number > 0) {
      curr_zi_order = feed_zi_ord;
    }

    if (!done) {
      std::vector<uint> tmp_vec;

      if (n > 1)
        tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end());

      // Compare if the food is the expected food; if not, store it for later use
      if (feed_zi_ord == tmp_vec) {
        // first check if we are done. If not start the reconstruction again using
        // the chinese remainder theorem in combining the previous results
        if (new_prime) {
          ti.emplace_back(new_ti);
          sub_num.clear();
          sub_den.clear();

          if (rec_rat_coef()) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              done = test_guess(num);
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

          new_prime = false;
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

            //todo remove if check = true
            if (!check)
              ai.emplace_back(comp_ai(i, i, num));
          }
        } else {
          if (coef_mat.empty())
            coef_mat.reserve(num_eqn);

          std::vector<std::pair<FFInt, FFInt>> t_food = {std::make_pair(new_ti, num)};

          // Prepare food for Gauss system
          if (n > 1) {
            std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end());

            if (saved_ti.find(tmp_vec) != saved_ti.end()) {
              t_food.insert(t_food.end(), saved_ti.at(tmp_vec).begin(), saved_ti.at(tmp_vec).end());
              saved_ti.erase(tmp_vec);
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
              }
            }

            if (prime_number == 0)
              yis.insert(yis.begin(), FFInt(1));
            else
              yis.insert(yis.begin(), tmp_ti);

            // build Gauss system for univariate reconstruction needed for
            // multivariate rational functions
            if (prime_number == 0)
              build_uni_gauss(tmp_ti, tmp_num, yis);
            else
              build_multi_gauss(tmp_num, yis);

            if (coef_mat.size() == num_eqn) {
              check = true;
              break;
            }
          }
        }

        if (check) {
          check = false;

          std::pair<PolynomialFF, PolynomialFF> canonical;

          // If the maximal/minimal degree of the polynomials are not set
          // determine them and save all information
          if (max_deg_num == -1) {
            if (ai.capacity() != ai.size()) {
              ai.shrink_to_fit();
              ti.shrink_to_fit();
            }

            ti.pop_back();

            canonical = construct_canonical();
            PolynomialFF denominator = canonical.second;

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

            max_deg_num = canonical.first.max_deg()[0];
            max_deg_den = canonical.second.max_deg()[0];

            if (n == 1) min_deg_den_vec = {canonical.second.max_deg()[0]};

            curr_deg_num = max_deg_num;

            if (max_deg_den > 0)
              curr_deg_den = max_deg_den;

            FFInt equializer = FFInt(1) / denominator.coefs[denominator.min_deg()];

            canonical.first = canonical.first * equializer;
            canonical.second = denominator * equializer;

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
            canonical = solve_multi_gauss();

          if (n == 1) {
            std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(canonical);
            combine_primes(tmp);
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              prime_number++;
              queue.clear();
            }
            saved_ti.clear();
            new_prime = true;
            return;
          } else if (prime_number == 0) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi = curr_zi;
            }

            ff_map num_coef = canonical.first.coefs;
            ff_map den_coef = canonical.second.coefs;

            // save the current results to the map to access them later
            for (int i = 0; i <= (int)(max_deg_num - tmp_solved_coefs_num); i++) {
              if (first_run) {
                {
                  std::unique_lock<std::mutex> lock_statics(mutex_statics);
                  PolyReconst rec(n - 1, (uint) i, true);
                  rec.set_anchor_points(anchor_points);
                  coef_n.emplace(std::make_pair((uint) i, std::move(rec)));
                }

                if (i < max_deg_num) {
                  std::vector<uint> zero_deg(n);
                  ff_map zero_mon = {{zero_deg, 0}};
                  sub_num.emplace(std::make_pair(i, PolynomialFF(n, zero_mon)));
                }
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
                  rec.set_anchor_points(anchor_points);
                  coef_d.emplace(std::make_pair(i, std::move(rec)));
                }

                if ((int) i < max_deg_den) {
                  std::vector<uint> zero_deg(n);
                  ff_map zero_mon = {{zero_deg, 0}};
                  sub_den.emplace(std::make_pair(i, PolynomialFF(n, zero_mon)));

                  if (i == 1)
                    sub_den.emplace(std::make_pair(0, PolynomialFF(n, zero_mon)));
                }
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

              for (const auto & el : denominator.coefs) {
                std::vector<uint> degs_reverse = el.first;
                std::reverse(degs_reverse.begin(), degs_reverse.end());

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

              // if the denominator is just a constant, there is no corresponding
              // PolyReconst object. Thus we set the minimal degree to a zero tuple
              FFInt const_shift = sub_den[0].calc(std::vector<FFInt> (n, 0));

              if (min_deg_den_vec.empty() || const_shift != 1) {
                min_deg_den_vec = std::vector<uint> (n);
                ff_map dummy_map;
                dummy_map.emplace(std::make_pair(min_deg_den_vec, FFInt(1) - const_shift));
                denominator = denominator + PolynomialFF(n, dummy_map);
              }

              coef_n.clear();
              coef_d.clear();
              curr_zi_order_num.clear();
              curr_zi_order_den.clear();

              FFInt first_coef = denominator.coefs[min_deg_den_vec];

              // normalize
              FFInt equializer = FFInt(1) / first_coef;

              numerator = numerator * equializer;
              denominator = denominator * equializer;

              std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(std::make_pair(numerator, denominator));

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
              for (uint zi = 2; zi <= n; zi ++) {
                auto key = std::make_pair(zi, curr_zi_order[zi - 2]);
                std::unique_lock<std::mutex> lock_statics(mutex_statics);
                set_new_rand(lock_statics, key);
              }
            }

            /*
             * If not finished, check if we can use some saved runs
             */
            std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end());

            if (saved_ti.find(tmp_vec) != saved_ti.end()) {
              std::pair<FFInt, FFInt> key_val = saved_ti.at(tmp_vec).back();
              saved_ti.at(tmp_vec).pop_back();
              interpolate(key_val.first, key_val.second, tmp_vec);
            }

            return;
          } else {
            if (n > 1) {
              std::unique_lock<std::mutex> lock(mutex_status);
              std::fill(curr_zi_order.begin(), curr_zi_order.end(), 1);
            }

            std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(canonical);
            combine_primes(tmp);
            std::unique_lock<std::mutex> lock(mutex_status);
            prime_number ++;
            new_prime = true;
          }
        }
      } else if (n > 1 && prime_number == 0 && feed_zi_ord != tmp_vec) {
        if (saved_ti.find(feed_zi_ord) == saved_ti.end()) {
          std::vector<std::pair<FFInt, FFInt>> tmp_ti = {std::make_pair(new_ti, num)};
          saved_ti[feed_zi_ord] = tmp_ti;
        } else {
          saved_ti.at(feed_zi_ord).emplace_back(std::make_pair(new_ti, num));
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
      //if (clock_test == 0) clock_test = clock();

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
            //std::cout << " adding shift\n";
            //std::clock_t begin = std::clock();
            PolynomialFF sub_pol = rec.get_result_ff().homogenize(curr_deg).add_shift(shift);
            sub_pol -= rec.get_result_ff().homogenize(curr_deg);
            //std::cout << "time : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";

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

        //std::cout << "deg " << curr_deg << " in num " << is_num << " took : " << float(clock() - clock_test) / CLOCKS_PER_SEC << "\n";

        /*
         * Remove already solved coefficients from Gauss eliminiation
         */
        curr_deg--;

        //clock_test = 0;

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
            non_solved_degs_num.emplace_back(c_ni.first);
        } catch (std::exception& e) {
          non_solved_degs_num.emplace_back(c_ni.first);
        }
      }

      mpz_map combined_di_back = combined_di;

      for (auto & c_di : combined_di_back) {
        try {
          RationalNumber rn = get_rational_coef(c_di.second, combined_prime);

          if (rn.numerator == c_di.second && rn.denominator == 1)
            remove_di(c_di.first, rn);
          else
            non_solved_degs_den.emplace_back(c_di.first);
        } catch (std::exception& e) {
          non_solved_degs_den.emplace_back(c_di.first);
        }
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

        try {
          RationalNumber last_rn = get_rational_coef(c_ni.second, combined_prime_back);
          RationalNumber curr_rn = get_rational_coef(combined_ni[c_ni.first], combined_prime);

          if (last_rn == curr_rn)
            remove_ni(c_ni.first, curr_rn);
          else
            non_solved_degs_num.emplace_back(c_ni.first);
        } catch (std::exception& e) {
          //TODO how does this work?
          if (c_ni.second == combined_ni[c_ni.first]) {
            RationalNumber rn = RationalNumber(c_ni.second, 1);
            remove_ni(c_ni.first, rn);
          } else
            non_solved_degs_num.emplace_back(c_ni.first);
        }
      }

      for (auto & c_di : combined_di_back) {
        try {
          RationalNumber last_rn = get_rational_coef(c_di.second, combined_prime_back);
          RationalNumber curr_rn = get_rational_coef(combined_di[c_di.first], combined_prime);

          if (last_rn == curr_rn)
            remove_di(c_di.first, curr_rn);
          else
            non_solved_degs_den.emplace_back(c_di.first);
        } catch (std::exception& e) {
          //TODO how does this work?

          if (c_di.second == combined_di[c_di.first]) {
            RationalNumber rn = RationalNumber(c_di.second, 1);
            remove_di(c_di.first, rn);
          } else
            non_solved_degs_den.emplace_back(c_di.first);
        }
      }

      combined_ni_back.clear();
      combined_di_back.clear();
      combined_prime_back = 0;
    }

    // Sort non solved coefficients to have a uniquely defined system of equations
    std::sort(non_solved_degs_num.begin(), non_solved_degs_num.end());
    std::sort(non_solved_degs_den.begin(), non_solved_degs_den.end());
    std::unique_lock<std::mutex> lock(mutex_status);
    num_eqn = non_solved_degs_num.size() + non_solved_degs_den.size();
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
      throw std::runtime_error("Access to unfinished result");
    }
  }

  bool RatReconst::rec_rat_coef() {
    bool run_test = true;

    for (const auto & ci : combined_ni) {
      mpz_class a = ci.second;

      try {
        g_ni[ci.first] = get_rational_coef(a, combined_prime);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    for (const auto & ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di[ci.first] = get_rational_coef(a, combined_prime);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    if (!run_test) {
      for (const auto & ci : combined_ni) {
        g_ni.erase(ci.first);
      }

      for (const auto & ci : combined_di) {
        g_di.erase(ci.first);
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

  std::pair<PolynomialFF, PolynomialFF> RatReconst::solve_gauss() {
    //std::clock_t begin = clock();
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

    //for(auto iii : results) std::cout << iii << "\n";
    //std::cout << " univ gauss (n = " << num_eqn << ") time : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";
    return std::make_pair(PolynomialFF(1, numerator), PolynomialFF(1, denominator));
  }

  std::pair< PolynomialFF, PolynomialFF > RatReconst::solve_multi_gauss() {
    std::vector<FFInt> results = solve_gauss_system(num_eqn, coef_mat);
    coef_mat.clear();

    // Bring result in canonical form
    ff_map numerator;
    ff_map denominator;

    int terms_num = non_solved_degs_num.size();

    for (int i = 0; i < terms_num; i ++) {
      std::vector<uint> power = non_solved_degs_num[i];
      numerator.emplace(std::make_pair(std::move(power), results[i]));
    }

    uint terms_den = non_solved_degs_den.size();

    for (uint i = 0; i < terms_den; i ++) {
      std::vector<uint> power = non_solved_degs_den[i];
      denominator.emplace(std::make_pair(std::move(power), results[i + terms_num]));
    }

    return std::make_pair(PolynomialFF(n, numerator), PolynomialFF(n, denominator));
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

    for (const auto coef : rf.first.coefs) {
      ci_mpz_1.emplace(coef.first, mpz_class(coef.second.n));
    }

    for (const auto coef : rf.second.coefs) {
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
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yis[i] = rand_zi[std::make_pair(i + 1, curr_zi_order[i - 1])];
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
    }

    return *this;
  }

  void RatReconst::disable_shift() {
    shift = std::vector<FFInt> (n, 0);
  }

  void RatReconst::build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);
    std::vector<FFInt> solved_coef_sub_num {};
    std::vector<FFInt> solved_coef_sub_den {};
    std::vector<FFInt> yis_wo_t = yis;
    yis_wo_t.erase(yis_wo_t.begin());

    // if(clock_test_2 == 0) clock_test_2 = clock();
    for (int r = 0; r <= max_deg_num; r++) {
      // If the current degree is smaller than the total degree of the polynomial
      // subtract the higher terms to save numerical runs
      if (shift[0] != 0 && coef_n.size() > 0 && coef_n[r].is_new_prime()) {
        FFInt sub;

        if (r < max_deg_num)
          sub = (coef_n[r].get_result_ff().calc(yis_wo_t) + sub_num[r].calc(yis)) * tmp_ti.pow(r);
        else
          sub = (coef_n[r].get_result_ff().calc(yis_wo_t)) * tmp_ti.pow(r);

        solved_coef_sub_num.emplace_back(sub);
      } else
        eq.emplace_back(tmp_ti.pow(FFInt(r)));
    }

    for (int rp = 1; rp <= max_deg_den; rp++) {
      // If the current degree is smaller than the total degree of the polynomial
      // subtract the higher terms to save numerical runs
      if (shift[0] != 0  && coef_d.size() > 0 && coef_d[rp].is_new_prime()) {
        FFInt sub;

        if (rp < max_deg_den)
          sub = (coef_d[rp].get_result_ff().calc(yis_wo_t) + sub_den[rp].calc(yis)) * tmp_ti.pow(rp);
        else
          sub = (coef_d[rp].get_result_ff().calc(yis_wo_t)) * tmp_ti.pow(rp);

        solved_coef_sub_den.emplace_back(sub);
      } else
        eq.emplace_back(-tmp_ti.pow(rp) * tmp_num);
    }

    // The lowest degree in univariate Gauss for multivariate polynomial
    // reconstruction is always zero. Hence, we just need to emplace the
    // evaluation of f(yis) at this point
    eq.emplace_back(tmp_num);

    for (auto & solved_coef_num : solved_coef_sub_num) {
      eq.back() += -solved_coef_num;
    }

    for (auto & solved_coef_den : solved_coef_sub_den) {
      eq.back() += solved_coef_den * tmp_num;
    }

    //std::cout << "time uni gauss : " << float(clock() - clock_test_2) / CLOCKS_PER_SEC << "\n";

    //clock_test_2 = 0;
    coef_mat.emplace_back(std::move(eq));
  }

  void RatReconst::build_multi_gauss(const FFInt& tmp_num, const std::vector<FFInt>& yis) {
    std::vector<FFInt> eq;
    eq.reserve(num_eqn + 1);

    // Increase the whole zi_order by 1
    if (n > 1) {
      {
        std::unique_lock<std::mutex> lock(mutex_status);
        std::transform(curr_zi_order.begin(), curr_zi_order.end(),
        curr_zi_order.begin(), [](uint x) {return x + 1;});
      }

    }

    // Build system of equations; in combined_.. are the non-solved coefficients
    for (const auto & pow_vec : non_solved_degs_num) {
      FFInt coef = 1;

      for (uint i = 0; i < n; i++) {
        coef *= yis[i].pow(pow_vec[i]);
      }

      eq.emplace_back(coef);
    }

    for (const auto & pow_vec : non_solved_degs_den) {
      FFInt coef = FFInt(0) - tmp_num;

      for (uint i = 0; i < n; i++) {
        coef *= yis[i].pow(pow_vec[i]);
      }

      eq.emplace_back(coef);
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

    for (const auto & el : g_di) {
      mpz_class tmp_1 = el.second.numerator % FFInt::p;
      mpz_class tmp_2 = el.second.denominator % FFInt::p;

      if (tmp_1 < 0) tmp_1 = FFInt::p + tmp_1;

      FFInt coef = FFInt(std::stoull(tmp_1.get_str())) / FFInt(std::stoull(tmp_2.get_str()));

      for (uint i = 0; i < n; i++) {

        coef *= yis[i].pow(el.first[i]);
      }

      eq.back() += coef * tmp_num;
    }

    coef_mat.emplace_back(std::move(eq));
  }

  uint64_t RatReconst::find_sieve_size(uint n) {
    // For small n, the formula returns a value too low, so we can just
    // hardcode the sieve size to 5 (5th prime is 11).
    if (n < 6)
      return 13;

    // We can't find a prime that will exceed ~0UL.
    if (n >= (~0UL / std::log(~0UL)))
      return 0;

    // Binary search for the right value.
    unsigned long low  = n;
    unsigned long high = ~0UL - 1;

    do {
      unsigned long mid   = low + (high - low) / 2;
      double        guess = mid / std::log(mid);

      if (guess > n)
        high = (unsigned long) mid - 1;
      else
        low = (unsigned long) mid + 1;
    } while (low < high);

    return high + 1;
  }

  uint64_t RatReconst::find_nth_prime(uint n) {
    if (!n) return 1;           // "0th prime"

    if (!--n) return 2;         // first prime

    unsigned long sieve_size = find_sieve_size(n);
    unsigned long count     = 0;
    unsigned long max_i     = sqrt(sieve_size - 1) + 1;

    if (sieve_size == 0)
      return 0;

    std::vector<bool> sieve(sieve_size);

    for (unsigned long i = 3;  true;  i += 2) {
      if (!sieve[i]) {
        if (++count == n)
          return i;

        if (i >= max_i)
          continue;

        unsigned long j    = i * i;
        unsigned long inc  = i + i;
        unsigned long maxj = sieve_size - inc;

        // This loop checks j before adding inc so that we can stop
        // before j overflows.
        do {
          sieve[j] = true;

          if (j >= maxj)
            break;

          j += inc;
        } while (1);
      }
    }

    return 0;
  }

  void RatReconst::generate_anchor_points(uint max_order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    gen_anchor_points(lock_statics, max_order);
  }
}
