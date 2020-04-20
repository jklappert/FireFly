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

#include "firefly/PolyReconst.hpp"
#include "firefly/Logger.hpp"
#include "firefly/Poly.hpp"
#include "firefly/ReconstHelper.hpp"
#include "firefly/utils.hpp"

namespace firefly {
  // TODO for new prime just use Vandermonde matrices to solve interpolation problem

  std::mutex PolyReconst::mutex_anchor;
  std::vector<FFInt> PolyReconst::global_anchor_points;

  PolyReconst::PolyReconst(uint32_t n_, const int deg_inp, const bool with_rat_reconst_inp) {
    std::lock_guard<std::mutex> lock_status(mutex_status);

    type = POLY;
    n = n_;
    zero_element = std::vector<uint32_t>(n);
    combined_prime = FFInt::p;
    curr_zi_order = std::vector<uint32_t>(n, 1);

    deg = deg_inp;
    with_rat_reconst = with_rat_reconst_inp;

    for (uint32_t i = 1; i <= n; ++i) {
      ais.emplace(std::make_pair(std::vector<uint32_t> (n), std::vector<FFInt> ()));
      max_deg.emplace(std::make_pair(i, -1));
    }

    if (!is_rand_zi_empty()) {
      private_anchor_points = global_anchor_points;
    }
  }

  PolyReconst::PolyReconst() {}

  void PolyReconst::set_anchor_points(const std::vector<FFInt>& anchor_points, bool force) {
    std::lock_guard<std::mutex> lock_statics(mutex_anchor);

    if (global_anchor_points.empty() || force) {
      global_anchor_points = std::vector<FFInt> (n, 0);

      for (uint32_t i = 0; i != n; ++i) {
	global_anchor_points[i] = anchor_points[i];
      }

      private_anchor_points = global_anchor_points;
    }
  }

  void PolyReconst::feed(const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t fed_prime) {
    std::lock_guard<std::mutex> lock(mutex_status);

    if (fed_prime == prime_number)
      queue.emplace_back(std::make_tuple(num, feed_zi_ord));
  }

  // Only use together with RatReconst, please.
  void PolyReconst::feed(const FFInt& num) {
    interpolate(num, curr_zi_order);
  }

  void PolyReconst::interpolate() {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (is_interpolating || queue.empty()) return;
    else {
      is_interpolating = true;

      while (!queue.empty()) {
        auto food = queue.front();
        queue.pop_front();
        lock.unlock();
        interpolate(std::get<0>(food), std::get<1>(food));
        lock.lock();
      }
    }

    is_interpolating = false;
  }

  void PolyReconst::interpolate(const FFInt& num, const std::vector<uint32_t>& feed_zi_ord) {
    if (feed_zi_ord == curr_zi_order) {
      if (!done) {
        if (new_prime && !with_rat_reconst) {
          // be sure that you have called generate_anchor_points!
          bool runtest = true;

          solved_degs.clear();

          ais.clear();
          nums_for_bt.clear();
          bt_terminator.clear();
          lambda.clear();
          b.clear();
          l.clear();
          delta.clear();
          bm_iteration.clear();

          ais.emplace(std::make_pair(std::vector<uint32_t> (n), std::vector<FFInt> ()));

          for (const auto ci : combined_ci) {
            mpz_class a = ci.second;

            auto res = get_rational_coef(a, combined_prime);

            if (res.first)
              gi.emplace(std::make_pair(ci.first, res.second));
            else {
              runtest = false;
              break;
            }
          }

          if (runtest) {
            {
              std::lock_guard<std::mutex> lock(mutex_status);
              done = test_guess(num);
            }

            if (done) {
              combined_prime = 0;
              combined_ci.clear();
              max_deg.clear();
              use_chinese_remainder = false;
              std::lock_guard<std::mutex> lock(mutex_status);
              new_prime = false;
              return;
            }
          }

          gi.clear();

          if (!use_chinese_remainder) use_chinese_remainder = true;

          std::lock_guard<std::mutex> lock(mutex_status);
          zi = 1;
          new_prime = false;
        }

        uint32_t i = curr_zi_order[zi - 1] - 1;

        // Univariate Newton interpolation for the lowest stage.
        if (zi == 1) {

          bool finished = false;

          if (use_newton) {
            if (i == 0)
              ais[zero_element].emplace_back(num);
            else
              ais[zero_element].emplace_back(comp_ai(i, num, ais[zero_element]));

            if (ais[zero_element].back() == 0) {
              combine_res = true;
              ais[zero_element].pop_back();
              finished = true;
            } else if (deg != -1 && static_cast<uint32_t>(deg) == i) {

              combine_res = true;
              finished = true;
            } else if (is_set_individual_degree_bounds == true && individual_degree_bounds[zi - 1] == i) {
              combine_res = true;
              finished = true;
	    }
          }

          if (use_bt && !finished) {
            if (nums_for_bt[zero_element].size() == 0 || i == 0) {
              nums_for_bt.erase(zero_element);
              lambda[zero_element] = {1};
              bt_terminator[zero_element] = 1;
              b[zero_element] = {0};
              l[zero_element] = 0;
              delta[zero_element] = 1;
              bm_iteration[zero_element] = 1;
            }

            nums_for_bt[zero_element].emplace_back(num);

            finished = berlekamp_massey_step(zero_element);

            if (finished) {
              std::pair<std::vector<FFInt>, std::vector<size_t>> roots = rootsexponents(zero_element, get_rand_zi(zi, 1));

              if (roots.first.size() == lambda[zero_element].size() - 1) {
                std::vector<FFInt> result;
                result = solve_transposed_vandermonde(roots.first, nums_for_bt[zero_element]);
                //combine the result if suceeded with the former interpolated polynomial
                //check for tmp solved degrees
                check_for_tmp_solved_degs_for_bt(zero_element, result, roots.second);
                combine_res = true;
                ais.erase(zero_element);
                bt_terminator.erase(zero_element);
                b.erase(zero_element);
                l.erase(zero_element);
                delta.erase(zero_element);
                bm_iteration.erase(zero_element);
                nums_for_bt.erase(zero_element);
                lambda.erase(zero_element);
              } else
                bt_terminator[zero_element] = 1;
            }
          } else if (use_bt && finished) {
            bt_terminator.erase(zero_element);
            b.erase(zero_element);
            l.erase(zero_element);
            delta.erase(zero_element);
            bm_iteration.erase(zero_element);
            nums_for_bt.erase(zero_element);
            lambda.erase(zero_element);
          }

          std::lock_guard<std::mutex> lock(mutex_status);
          curr_zi_order[zi - 1]++;
        } else {
          // Build Vandermonde system
          FFInt res = num;

          for (const auto & el : solved_degs) {
            std::vector<uint32_t> deg_vec = el.first;
            FFInt coef_num = el.second;

            for (uint32_t tmp_zi = 1; tmp_zi <= zi; ++tmp_zi) {
              // curr_zi_ord starts at 1, thus we need to subtract 1 entry
              coef_num *= get_rand_zi(tmp_zi, curr_zi_order[tmp_zi - 1]).pow(deg_vec[tmp_zi - 1]);
            }

            res -= coef_num;
          }

          for (const auto & el : tmp_solved_degs) {
            std::vector<uint32_t> deg_vec = el.first;
            FFInt coef_num = el.second;

            for (uint32_t tmp_zi = 1; tmp_zi <= zi; ++tmp_zi) {
              // curr_zi_ord starts at 1, thus we need to subtract 1 entry
              coef_num *= get_rand_zi(tmp_zi, curr_zi_order[tmp_zi - 1]).pow(deg_vec[tmp_zi - 1]);
            }

            res -= coef_num;
          }

          nums.emplace_back(res);

          // to total degree (save them in solved degs including their coefficient
          // to subtract them)
          // Solve Vandermonde system and calculate the next a_i
          if (nums.size() == rec_degs.size()) {
            const uint32_t order_save = curr_zi_order[zi - 1];
            {
              std::lock_guard<std::mutex> lock(mutex_status);
              curr_zi_order = std::vector<uint32_t> (n, 1);
              curr_zi_order[zi - 1] = order_save + 1;
            }

            uint32_t not_done_counter_newton = 0;
            uint32_t not_done_counter_bt = 0;

            for (const auto & el : build_and_solve_transposed_vandermonde()) {
              std::vector<uint32_t> key = el.first;

              bool finished = false;

              if (use_newton) {
                uint32_t tmp_deg = deg;

                for (auto ele : key) {
                  tmp_deg -= ele;
                }

                FFInt tmp_ai = comp_ai(i, el.second, ais[key]);
                ais[key].emplace_back(tmp_ai);

                if (tmp_ai == 0) {
                  ais[key].pop_back();
                  finished = true;
                } else if (deg != -1 && i == tmp_deg) {
                  finished = true;
                }  else if (is_set_individual_degree_bounds == true && individual_degree_bounds[zi - 1] == i) {
		  finished = true;
		} else
                  ++not_done_counter_newton;

                if (finished) {
                  check_for_tmp_solved_degs_for_newton(key, ais[key]);
                  ais.erase(key);

                  if (use_bt) {
                    bt_terminator.erase(key);
                    b.erase(key);
                    l.erase(key);
                    delta.erase(key);
                    bm_iteration.erase(key);
                    nums_for_bt.erase(key);
                    lambda.erase(key);
                  }
                }
              }

              if (use_bt && !finished) {
                nums_for_bt[key].emplace_back(el.second);

                if (nums_for_bt[key].size() == 1) {
                  lambda[key] = {1};
                  bt_terminator[key] = 1;
                  b[key] = {0};
                  l[key] = 0;
                  delta[key] = 1;
                  bm_iteration[key] = 1;
                }

                finished = berlekamp_massey_step(key);

                if (finished) {
                  std::pair<std::vector<FFInt>, std::vector<size_t>> roots = rootsexponents(key, get_rand_zi(zi, 1));

                  if (roots.first.size() == lambda[key].size() - 1) {
                    std::vector<FFInt> result = solve_transposed_vandermonde(roots.first, nums_for_bt[key]);
                    //combine the result if suceeded with the former interpolated polynomial
                    //check for tmp solved degrees
                    check_for_tmp_solved_degs_for_bt(key, result, roots.second);
                    bt_terminator.erase(key);
                    b.erase(key);
                    l.erase(key);
                    delta.erase(key);
                    bm_iteration.erase(key);
                    nums_for_bt.erase(key);
                    lambda.erase(key);
                    ais.erase(key);
                  } else {
                    ++not_done_counter_bt;
                    bt_terminator[key] = 1;
                  }
                } else
                  ++not_done_counter_bt;
              }
            }

            if (not_done_counter_newton == 0 && use_newton)
              combine_res = true;

            if (not_done_counter_bt == 0 && use_bt && !combine_res)
              combine_res = true;

          } else {
            // increase all zi order of the lower stages by one
            std::lock_guard<std::mutex> lock(mutex_status);

            for (uint32_t tmp_zi = 1; tmp_zi < zi; ++tmp_zi) {
              curr_zi_order[tmp_zi - 1]++;
            }
          }
        }

        // if the lowest stage ai is zero, combine them into an ai for a higher stage
        // and check if we are done
        if (combine_res) {
          combine_res = false;

          if (n > 1) {
            // combine the current stage with the multivariate polynomial of the
            // previous stages and extract the reconstructed degrees to prepare
            // the gauss system
            // Remove all terms which are of total degree of the polynomial
            // to remove them from the next Vandermonde systems
            rec_degs.clear();

            if (zi == 1 && use_newton)
              check_for_tmp_solved_degs_for_newton(zero_element, ais[zero_element]);

            ff_map pol_ff = tmp_solved_degs;

            tmp_solved_degs.clear();

            for (auto & el : pol_ff) {
              rec_degs.emplace_back(el.first);
            }

            if (rec_degs.size() == 0 && zi != n) {
              std::lock_guard<std::mutex> lock(mutex_status);
              zi = n;
            }

            if (zi != n) {
              std::lock_guard<std::mutex> lock(mutex_status);
              zi ++;

              nums.reserve(rec_degs.size());
              ais.clear();

              nums_for_bt.clear();
              bt_terminator.clear();
              b.clear();
              l.clear();
              lambda.clear();
              delta.clear();
              bm_iteration.clear();

              for (const auto & el : pol_ff) {
                ais[el.first] = {el.second};
                nums_for_bt[el.first] = {el.second};
                lambda[el.first] = {1};
                bt_terminator[el.first] = 1;
                b[el.first] = {0};
                l[el.first] = 0;
                delta[el.first] = 1;
                bm_iteration[el.first] = 1;
                berlekamp_massey_step(el.first);
              }

              // reset zi order
              curr_zi_order = std::vector<uint32_t> (n, 1);
              curr_zi_order[zi - 1] = 2;
            } else {
              check = true;
              ais.clear();
              nums_for_bt.clear();
              bt_terminator.clear();
              lambda.clear();
              b.clear();
              l.clear();
              delta.clear();
              bm_iteration.clear();

              for (const auto & el : pol_ff) {
                ais[el.first].emplace_back(el.second);
                nums_for_bt[el.first].emplace_back(el.second);
                lambda[el.first] = {1};
                bt_terminator[el.first] = 1;
                b[el.first] = {0};
                l[el.first] = 0;
                delta[el.first] = 1;
                bm_iteration[el.first] = 1;
                berlekamp_massey_step(el.first);
              }
            }
          } else if (zi == 1 && n == 1)
            check = true;

          if (check && zi == n) {
            {
              std::lock_guard<std::mutex> lock(mutex_status);
              curr_zi_order = std::vector<uint32_t> (n, 1);
            }

            ff_map tmp_pol_ff {};

            for (const auto & el : ais) {
              ff_map tmp = construct_tmp_canonical(el.first, el.second);
              tmp_pol_ff.insert(tmp.begin(), tmp.end());
            }

            tmp_pol_ff.insert(solved_degs.begin(), solved_degs.end());

            if (!with_rat_reconst) {
              mpz_map ci_tmp = convert_to_mpz(tmp_pol_ff);

              if (!use_chinese_remainder) {
                combined_ci = ci_tmp;
              } else {
                // use another prime to utilize the Chinese Remainder Theorem to reconstruct the rational
                // coefficients

                std::pair<mpz_class, mpz_class> p1;
                std::pair<mpz_class, mpz_class> p2;
                std::pair<mpz_class, mpz_class> p3;

                for (const auto el : combined_ci) {
                  if (ci_tmp.find(el.first) == ci_tmp.end() && gi.find(el.first) == gi.end()) {
                    ci_tmp.emplace(std::make_pair(el.first, 0));
                  }
                }

                for (auto it = ci_tmp.begin(); it != ci_tmp.end(); ++it) {
                  p2 = std::make_pair(it->second, FFInt::p);

                  if (combined_ci.find(it->first) == combined_ci.end() && gi.find(it->first) == gi.end()) {
                    combined_ci.emplace(std::make_pair(it->first, 0));
                  }

                  p1 = std::make_pair(combined_ci[it->first], combined_prime);
                  p3 = run_chinese_remainder(p1, p2);
                  combined_ci[it->first] = p3.first;
                }

                combined_prime = p3.second;
              }
            } else {
              ais.clear();

              nums_for_bt.clear();
              bt_terminator.clear();
              lambda.clear();
              b.clear();
              l.clear();
              delta.clear();
              bm_iteration.clear();

              result_ff = PolynomialFF(n, tmp_pol_ff).homogenize(deg);
              result_ff.n = n + 1;
              rec_degs = std::vector<std::vector<uint32_t>>();
              nums = std::vector<FFInt>();
            }

            std::lock_guard<std::mutex> lock(mutex_status);
            new_prime = true;
            prime_number ++;
            check = false;
            return;
          }
        }
      }
    }
  }

  Polynomial PolyReconst::get_result() {
    if (result.coefs.empty()) {
      if (done) {
        result = Polynomial(gi);
        result.sort();
        gi.clear();
      }
    }

    return result;
  }

  PolynomialFF PolyReconst::get_result_ff() {
    return result_ff;
  }

  FFInt PolyReconst::comp_ai(int i, const FFInt& num, const std::vector<FFInt>& ai) const {
    FFInt res = num;

    if (i > 0) {
      FFInt yi_i_p_1 = get_rand_zi(zi, i + 1);

      for (int ip_tmp = 1; ip_tmp != i + 1; ip_tmp++) {
        FFInt yi_ip = get_rand_zi(zi, ip_tmp);
        res = (res - ai[ip_tmp - 1]) / (yi_i_p_1 - yi_ip); //yi[i + 1] - yi[ip - 1 + 1] the +1 in the first and the +1 in the second is due to the 0th element in the vector*/
      }
    }

    return res;
  }

  ff_map PolyReconst::construct_canonical(const std::vector<FFInt>& ai) const {
    size_t size = ai.size();

    if (size == 1) return {{std::vector<uint32_t> (1, 0), ai[0]}};
    else if (size == 0) return {{std::vector<uint32_t> (n, 0), 0}};

    return (PolynomialFF(1, {{std::vector<uint32_t> (1, 0), ai[0]}}) + iterate_canonical(ai)).coefs;
  }

  PolynomialFF PolyReconst::iterate_canonical(const std::vector<FFInt>& ai) const {
    FFInt yi = get_rand_zi(zi, ai.size() - 1);

    PolynomialFF res = PolynomialFF(1, {{{0}, ai.back()}});
    res = res * (-yi) + res.mul(1);

    for (uint32_t i_tmp = ai.size() - 2; i_tmp != 0; i_tmp--) {
      yi = get_rand_zi(zi, i_tmp);
      res += PolynomialFF(1, {{{0}, ai[i_tmp]}});
      res = res * (-yi) + res.mul(1); //yi[i - 1 + 1] +1 is due to 0th element
    }

    return res;
  }

  bool PolyReconst::test_guess(const FFInt& num) {
    ff_map gi_ffi = convert_to_ffint(gi);
    PolynomialFF gy(n, gi_ffi);
    std::vector<FFInt> chosen_yi(n);

    for (uint32_t i = 1; i <= n; ++i) {
      chosen_yi[i - 1] = get_rand_zi(i, 1);
    }

    return gy.calc(chosen_yi) == num;
  }

  // Solves the Vandermonde linear system V*x=a
  // V is build from vis, x contain our coefficients, and a is the numerical
  // value of the function which should be interpolated for a given numerical
  // input
  ff_map PolyReconst::build_and_solve_transposed_vandermonde() {
    uint32_t num_eqn = rec_degs.size();
    std::vector<FFInt> result(num_eqn);

    // calculate base entries of Vandermonde matrix
    std::vector<FFInt> vis;
    vis.reserve(num_eqn);

    for (const auto & el : rec_degs) {
      FFInt vi = 1;

      for (uint32_t tmp_zi = 1; tmp_zi < zi; ++tmp_zi) {
        // curr_zi_ord starts at 1, thus we need to subtract 1 entry
        vi *= get_rand_zi(tmp_zi, el[tmp_zi - 1]);
      }

      vis.emplace_back(vi);
    }

    result = solve_transposed_vandermonde(vis, nums);

    // Bring result in canonical form
    ff_map poly;
    poly.reserve(num_eqn);

    for (uint32_t i = 0; i != num_eqn; ++i) {
      poly.emplace(std::make_pair(rec_degs[i], result[i]));
    }

    nums.clear();
    return poly;
  }

  void PolyReconst::generate_anchor_points() {
    std::lock_guard<std::mutex> lock(mutex_anchor);
    global_anchor_points = std::vector<FFInt> (n, 0);

    for (uint32_t i = 0; i != n; ++i) {
      global_anchor_points[i] = FFInt(get_rand_64());
    }

    private_anchor_points = global_anchor_points;
  }

  FFInt PolyReconst::get_rand_zi(uint32_t zi, uint32_t order) const {
    return private_anchor_points[zi - 1].pow(order);
  }

  std::vector<FFInt> PolyReconst::get_rand_zi_vec(const std::vector<uint32_t>& orders) const {
    std::vector<FFInt> yis {};

    for (uint32_t i = 0; i < n; ++i) {
      yis.emplace_back(private_anchor_points[i].pow(orders[i]));
    }

    return yis;
  }

  bool PolyReconst::is_rand_zi_empty() const {
    std::lock_guard<std::mutex> lock(mutex_anchor);
    return global_anchor_points.empty();
  }

  void PolyReconst::reset() {
    std::lock_guard<std::mutex> lock(mutex_anchor);
    global_anchor_points.clear();
    BaseReconst::reset();
  }

  ff_map PolyReconst::construct_tmp_canonical(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) const {
    ff_map tmp {};

    if (ai.size() == 1 && zi == n) {
      tmp.emplace(std::make_pair(deg_vec, ai[0]));
    } else {
      for (auto & el : construct_canonical(ai)) { // homogenize
        if (el.second != 0) {
          std::vector<uint32_t> new_deg(n);
          new_deg[zi - 1] = el.first[0];

          for (uint32_t j = 0; j < zi - 1; j++) {
            new_deg[j] = deg_vec[j];
          }

          tmp.emplace(std::make_pair(new_deg, el.second));
        }
      }
    }

    return tmp;
  }

  void PolyReconst::set_bt_threshold(size_t threshold) {
    bt_threshold = threshold;
  }

  bool PolyReconst::berlekamp_massey_step(const std::vector<uint32_t>& key) {
    FFInt delta_r = 0;
    size_t size_lambda = lambda[key].size();

    for (size_t i = 0; i != size_lambda; i++) {
      delta_r += lambda[key].at(i) * nums_for_bt[key].at(bm_iteration[key] - i - 1);
    }

    if (delta_r == 0) {
      b[key].insert(b[key].begin(), FFInt(0));

      while (b[key].back() == 0) {
        b[key].pop_back();
      }

      bt_terminator[key]++;

      if (bt_terminator[key] >= bt_threshold + 1 && bm_iteration[key] > 2 * l[key]) {
        bm_iteration[key]++;
        //return false;
        return true;
      } else {
        bm_iteration[key]++;
        return false;
      }
    }

    if (delta_r != 0) {
      std::vector<FFInt> b_temp(lambda[key]);
      std::vector<FFInt> lambda_temp;
      b[key].insert(b[key].begin(), 0);
      size_t size_b = b[key].size();

      for (size_t j = 0; j < size_lambda || j < size_b; ++j) {
        if (j < size_lambda && j < size_b)
          lambda_temp.emplace_back(lambda[key].at(j) - delta_r / delta[key]*b[key].at(j));

        if (j < size_lambda && j >= size_b)
          lambda_temp.emplace_back(lambda[key].at(j));

        if (j >= size_lambda && j < size_b)
          lambda_temp.emplace_back(- delta_r / delta[key]*b[key].at(j));
      }

      if (2 * l[key] < bm_iteration[key]) {
        l[key] = bm_iteration[key] - l[key];
        delta[key] = delta_r;
        b[key].swap(b_temp);
      }

      lambda[key].swap(lambda_temp);

      while (lambda[key].back() == 0) {
        lambda[key].pop_back();
      }

      while (b[key].back() == 0) {
        b[key].pop_back();
      }

      bm_iteration[key]++;

      bt_terminator[key] = 1;

      return false;
    }

    ERROR_MSG("Berlekamp/Massey step exit was wrong!");
    return false;
  }

  std::pair<std::vector<FFInt>, std::vector<size_t>> PolyReconst::rootsexponents(const std::vector<uint32_t>& key, const FFInt& base) {
    std::pair<std::vector<FFInt>, std::vector<size_t>> roots;
    FFInt a(1);
    size_t count = 0;
    size_t sol_deg = 0;

    for (size_t i = 0; i < key.size(); i++) {
      sol_deg += key.at(i);
    }

    while (roots.first.size() < lambda[key].size() - 1) {
      FFInt result(0);

      for (size_t i = 0; i < lambda[key].size(); i++) {
        result += lambda[key].at(lambda[key].size() - i - 1) * a.pow(i);
      }

      if (result == 0) {
        roots.first.emplace_back(a);
        roots.second.emplace_back(count);
      }

      a *= base;

      if (a == FFInt(1))
        break;

      count++;

      if (count > deg - sol_deg)
        break;
    }

    if (roots.first.size() != lambda[key].size() - 1)
      ERROR_MSG("The Polynomial calculated by the Berlekamp/Massey algorithm is not correct");

    return roots;
  }

  std::pair<std::vector<FFInt>, std::vector<size_t>> PolyReconst::rootsexponents_with_poly_class(const std::vector<uint32_t>& key, const FFInt& base) {
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponent;
    Poly lambdapoly(lambda[key]);
    lambdapoly.rev();
    std::vector<FFInt> roots = lambdapoly.roots();
    FFInt a(1);

    for (size_t count = 0; count <= (size_t) deg; count++) {
      for (size_t j = 0; j < roots.size(); j++) {
        if (a == roots.at(j)) {
          rootsexponent.first.emplace_back(a);
          rootsexponent.second.emplace_back(count);
        }
      }

      a *= base;

      if (a == FFInt(1))
        break;

      if (rootsexponent.first.size() == roots.size())
        break;
    }

    if (rootsexponent.first.size() != lambdapoly.get_deg())
      ERROR_MSG("The Polynomial calculated by the Berlekamp/Massey algorithm is not correct");

    return rootsexponent;
  }

  std::vector<FFInt> PolyReconst::solve_transposed_vandermonde(const std::vector<FFInt>& vis, const std::vector<FFInt>& fis) const {
    uint32_t num_eqn = vis.size();

    if (num_eqn == 0) {
      return std::vector<FFInt> {FFInt(0)};
    }

    std::vector<FFInt> result(num_eqn);
    // Initialize the coefficient vector of the master polynomial
    std::vector<FFInt> cis(num_eqn);

    // The coefficients of the master polynomial are found by recursion
    // where we have
    // P(Z) = (Z - v_0)*(Z - v_1)*...*(Z - v_{n-1})
    //      =  c_0 + c_1*Z + ... + Z^n
    cis[num_eqn - 1] = -vis[0];

    for (uint32_t i = 1; i < num_eqn; ++i) {
      for (uint32_t j = num_eqn - 1 - i; j < num_eqn - 1; ++j) {
        cis[j] -= vis[i] * cis[j + 1];
      }

      cis[num_eqn - 1] -= vis[i];
    }

    // Each subfactor in turn is synthetically divided,
    // matrix-multiplied by the right hand-side,
    // and supplied with a denominator (since all vi should be different,
    // there is no additional check if a coefficient in synthetical division
    // leads to a vanishing denominator)
    for (uint32_t i = 0; i < num_eqn; ++i) {
      FFInt t = 1;
      FFInt b = 1;
      FFInt s = fis[num_eqn - 1];

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += fis[j - 1] * b;
        t = vis[i] * t + b;
      }

      result[i] = s / t / vis[i];
    }

    return result;
  }

  void PolyReconst::check_for_tmp_solved_degs_for_newton(const std::vector<uint32_t>& deg_vec,
                                                         const std::vector<FFInt>& ai) {
    ff_map tmp = construct_tmp_canonical(deg_vec, ai);

    for (const auto & el : tmp) {
      int total_deg = 0;

      for (const auto & e : el.first) total_deg += e;

      if (total_deg == deg && el.second != 0)
        solved_degs.emplace(std::make_pair(el.first, el.second));
      else if (el.second != 0)
        tmp_solved_degs.emplace(std::make_pair(el.first, el.second));
    }

    if (zi > 1) {
      std::vector<std::vector<uint32_t>>::iterator it = std::find(rec_degs.begin(), rec_degs.end(), deg_vec);
      rec_degs.erase(it);
    }
  }

  void PolyReconst::check_for_tmp_solved_degs_for_bt(const std::vector<uint32_t>& deg_vec,
                                                     const std::vector<FFInt>& coeffs,
                                                     const std::vector<size_t>& exponents) {
    ff_map tmp;

    for (size_t i = 0; i != exponents.size(); ++i) {
      std::vector<uint32_t> new_deg_vec = deg_vec;
      new_deg_vec[zi - 1] = exponents[i];
      tmp.emplace(std::make_pair(new_deg_vec, coeffs[i]));
    }

    if (exponents.size() == 0)
      tmp.emplace(std::make_pair(deg_vec, 0));

    for (const auto & el : tmp) {
      int total_deg = 0;

      for (const auto & e : el.first) total_deg += e;

      if (total_deg == deg && el.second != 0)
        solved_degs.emplace(std::make_pair(el.first, el.second));
      else if (n == 1)
        solved_degs.emplace(std::make_pair(el.first, el.second));
      else if (el.second != 0)
        tmp_solved_degs.emplace(std::make_pair(el.first, el.second));
    }

    if (zi > 1) {
      std::vector<std::vector<uint32_t>>::iterator it = std::find(rec_degs.begin(), rec_degs.end(), deg_vec);
      rec_degs.erase(it);
    }
  }

  void PolyReconst::set_newton(bool use_newton_new) {
    use_newton = use_newton_new;
  }

  void PolyReconst::set_bt(bool use_bt_new) {
    use_bt = use_bt_new;
  }

  uint32_t PolyReconst::get_vandermonde_num_eqn() const {
    std::lock_guard<std::mutex> lock(mutex_status);

    if (zi == 1)
      return 1;
    else
      return rec_degs.size() - nums.size();
  }

  void PolyReconst::set_individual_degree_bounds(const std::vector<uint32_t>& individual_degree_bounds_) {
    is_set_individual_degree_bounds = true;
    individual_degree_bounds = individual_degree_bounds_;
  }
}
