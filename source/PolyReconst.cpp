// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================
#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"
#include <chrono>

namespace firefly {
  // TODO check if this interpolates in combination with RatReconst to use the
  // static rand_zi of RatReconst to save additional memory -> note that
  // zi order in RatReconst and PolyReconst is not the same!
  // TODO for new prime just use Vandermonde matrices to solve interpolation problem

  ff_pair_map PolyReconst::rand_zi;
  std::mutex PolyReconst::mutex_statics;

  PolyReconst::PolyReconst(uint32_t n_, const int deg_inp, const bool with_rat_reconst_inp) {
    std::unique_lock<std::mutex> lock_status(mutex_status);

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
  }

  PolyReconst::PolyReconst() {}

  void PolyReconst::set_anchor_points(const std::vector<FFInt>& anchor_points, bool force) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    if (rand_zi.empty() || force) {
      rand_zi.clear();

      for (uint32_t i = 1; i <= n; ++i) {
        rand_zi.emplace(std::make_pair(std::make_pair(i, 0), 1));
        rand_zi.emplace(std::make_pair(std::make_pair(i, 1), anchor_points[i - 1]));
      }
    }
  }

  void PolyReconst::feed(const FFInt& num, const std::vector<uint32_t>& feed_zi_ord, const uint32_t fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (fed_prime == prime_number)
      queue.emplace_back(std::make_tuple(num, feed_zi_ord));
  }

  // this function is not thread-safe; therefore, it should only be called from RatReconst
  void PolyReconst::feed(const std::vector<FFInt>& new_yis, const FFInt& num) {
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      for (uint32_t j = 0; j < n; ++j) {
        rand_zi.emplace(std::make_pair(j + 1, curr_zi_order[j]), new_yis[j]);
      }
    }

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
              std::unique_lock<std::mutex> lock(mutex_status);
              done = test_guess(num);
            }

            if (done) {
              combined_prime = 0;
              combined_ci.clear();
              max_deg.clear();
              use_chinese_remainder = false;
              std::unique_lock<std::mutex> lock(mutex_status);
              new_prime = false;
              return;
            }
          }

          gi.clear();

          if (!use_chinese_remainder) use_chinese_remainder = true;

          std::unique_lock<std::mutex> lock(mutex_status);
          zi = 1;
          new_prime = false;
        }

        uint32_t i = curr_zi_order[zi - 1] - 1;

        // Univariate Newton interpolation for the lowest stage.
        if (zi == 1) {
          if (i == 0)
            ais[zero_element].emplace_back(num);
          else {
            ais[zero_element].emplace_back(comp_ai(i, i, num, ais[zero_element]));

            if (ais[zero_element].back() == 0) {
              combine_res = true;
              ais[zero_element].pop_back();
            } else if (deg != -1 && (uint32_t) deg == i) {
              combine_res = true;
            }
          }

          std::unique_lock<std::mutex> lock(mutex_status);
          curr_zi_order[zi - 1] ++;
        } else {
          // Build Vandermonde system
          FFInt res = num;

          for (const auto & el : solved_degs) {
            std::vector<uint32_t> deg_vec = el.first;
            FFInt coef_num = el.second;

            for (uint32_t tmp_zi = 1; tmp_zi < zi; ++tmp_zi) {
              // curr_zi_ord starts at 1, thus we need to subtract 1 entry
              std::unique_lock<std::mutex> lock_statics(mutex_statics);
              coef_num *= rand_zi[std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 1])].pow(deg_vec[tmp_zi - 1]);
            }

            res -= coef_num;
          }

          for (const auto & el : tmp_solved_degs) {
            std::vector<uint32_t> deg_vec = el.first;
            FFInt coef_num = el.second;

            for (uint32_t tmp_zi = 1; tmp_zi < zi + 1; ++tmp_zi) {
              // curr_zi_ord starts at 1, thus we need to subtract 1 entry
              std::unique_lock<std::mutex> lock_statics(mutex_statics);
              coef_num *= rand_zi[std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 1])].pow(deg_vec[tmp_zi - 1]);
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
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = std::vector<uint32_t> (n, 1);
              curr_zi_order[zi - 1] = order_save + 1;
            }

            uint32_t not_done_counter = 0;

            for (const auto & el : solve_transposed_vandermonde()) {
              std::vector<uint32_t> key = el.first;

              FFInt tmp_ai = comp_ai(i, i, el.second, ais[key]);
              ais[key].emplace_back(tmp_ai);

              if (tmp_ai == 0) {
                ais[key].pop_back();
                check_for_tmp_solved_degs(key, ais[key]);
                ais.erase(key);
              } else if (deg != -1 && i == (uint32_t) deg) {
                check_for_tmp_solved_degs(key, ais[key]);
                ais.erase(key);
              } else {
                ++not_done_counter;
              }
            }

            if (not_done_counter == 0)
              combine_res = true;
          } else {
            // increase all zi order of the lower stages by one
            std::unique_lock<std::mutex> lock(mutex_status);

            for (uint32_t tmp_zi = 1; tmp_zi < zi; ++tmp_zi) {
              curr_zi_order[tmp_zi - 1] ++;
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

            if (zi == 1)
              check_for_tmp_solved_degs(zero_element, ais[zero_element]);

            ff_map pol_ff = tmp_solved_degs;
            tmp_solved_degs.clear();

            for (auto & el : pol_ff) {
              rec_degs.emplace_back(el.first);
              std::cout << "rec_degs " << deg << " " << el.first[0] << " " << el.first[1] << " "<< el.first[2] << "\n";
            }

            if (rec_degs.size() == 0 && zi != n) {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi = n;
            }

            if (zi != n) {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi ++;
              // The monomials which have to be reconstructed have to
              // ordered in a monotonical way to utilize the Vandermonde
              // system solver
              std::sort(rec_degs.begin(), rec_degs.end(), std::greater<std::vector<uint32_t>>());

              nums.reserve(rec_degs.size());
              ais.clear();

              for (const auto & el : pol_ff) {
                ais[el.first].emplace_back(el.second);
              }

              // reset zi order
              curr_zi_order = std::vector<uint32_t> (n, 1);
              curr_zi_order[zi - 1] = 2;
            } else {
              check = true;
              ais.clear();

              for (const auto & el : pol_ff) {
              std::cout << "fin_degs " << deg << " " << zi << " " << el.first[0] << " " << el.first[1] << " "<< el.first[2] << "\n";
                ais[el.first].emplace_back(el.second);
              }
            }
          } else if (zi == 1 && n == 1)
            check = true;

          if (check && zi == n) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = std::vector<uint32_t> (n, 1);
            }

            ff_map tmp_pol_ff {};

            for (const auto & el : ais) {
              std::cout << "1test " << deg << " " << el.first[0] << " " << el.first[1] << " "<< el.first[2] << "\n";;
              ff_map tmp = construct_tmp_canonical(el.first, el.second);
              std::cout << "test " << deg << " " << PolynomialFF(n, tmp);
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

                for (auto it = combined_ci.begin(); it != combined_ci.end(); ++it) {
                  p1 = std::make_pair(it->second, combined_prime);
                  p2 = std::make_pair(ci_tmp[it->first], FFInt::p);
                  p3 = run_chinese_remainder(p1, p2);
                  combined_ci[it->first] = p3.first;
                }

                combined_prime = p3.second;
              }
            } else {
              ais.clear();
              result_ff = PolynomialFF(n, tmp_pol_ff).homogenize(deg);
              result_ff.n = n + 1;
              std::cout << "res deg " << deg << " " << result_ff;
              rec_degs = std::vector<std::vector<uint32_t>>();
              nums = std::vector<FFInt>();
            }

            std::unique_lock<std::mutex> lock(mutex_status);
            new_prime = true;
            prime_number ++;
            check = false;
            return;
          }
        }
      }

      if (!with_rat_reconst) {
        for (uint32_t tmp_zi = 1; tmp_zi <= n; ++tmp_zi) {
          auto key = std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 1]);
          std::unique_lock<std::mutex> lock_statics(mutex_statics);

          if (rand_zi.find(key) == rand_zi.end())
            rand_zi.emplace(std::make_pair(key, rand_zi[std::make_pair(tmp_zi, 1)].pow(key.second)));
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

  FFInt PolyReconst::comp_ai(int i, int ip, const FFInt& num, std::vector<FFInt>& ai) {
    if (ip == 0) return num;

    FFInt yi_i_p_1;
    FFInt yi_ip;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yi_i_p_1 = rand_zi[std::make_pair(zi, i + 1)];
      yi_ip = rand_zi[std::make_pair(zi, ip)];
    }

    return (comp_ai(i, ip - 1, num, ai) - ai[ip - 1]) / (yi_i_p_1 - yi_ip); //yi[i + 1] - yi[ip - 1 + 1] the +1 in the first and the +1 in the second is due to the 0th element in the vector
  }

  ff_map PolyReconst::construct_canonical(const std::vector<FFInt>& ai) const {
    size_t size = ai.size();

    if (size == 1) return {{std::vector<uint32_t> (1, 0), ai[0]}};
    else if (size == 0) return {{std::vector<uint32_t> (n, 0), 0}};

    return (PolynomialFF(1, {{std::vector<uint32_t> (1, 0), ai[0]}}) + iterate_canonical(1, ai)).coefs;
  }

  PolynomialFF PolyReconst::iterate_canonical(uint32_t i, const std::vector<FFInt>& ai) const {
    FFInt yi;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yi = rand_zi[std::make_pair(zi, i)];
    }

    PolynomialFF dum_pol = PolynomialFF(1, {{{0}, ai[i]}});

    if (i < ai.size() - 1) {
      PolynomialFF poly = dum_pol + iterate_canonical(i + 1, ai);
      return poly.mul(1) + poly * (-yi); //yi[i - 1 + 1] +1 is due to 0th element
    }

    return dum_pol * (-yi) + dum_pol.mul(1); // yi[i - 1 + 1]
  }

  bool PolyReconst::test_guess(const FFInt& num) {
    ff_map gi_ffi = convert_to_ffint(gi);
    PolynomialFF gy(n, gi_ffi);
    std::vector<FFInt> chosen_yi(n);

    for (uint32_t i = 1; i <= n; ++i) {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      chosen_yi[i - 1] = rand_zi[std::make_pair(i, 1)];
    }

    return gy.calc(chosen_yi) == num;
  }

  // Solves the Vandermonde linear system V*x=a
  // V is build from vis, x contain our coefficients, and a is the numerical
  // value of the function which should be interpolated for a given numerical
  // input
  ff_map PolyReconst::solve_transposed_vandermonde() {
    uint32_t num_eqn = rec_degs.size();
    std::vector<FFInt> result(num_eqn);

    // calculate base entries of Vandermonde matrix
    std::vector<FFInt> vis;
    vis.reserve(num_eqn);

    for (const auto & el : rec_degs) {
      FFInt vi = 1;

      for (uint32_t tmp_zi = 1; tmp_zi < zi; ++tmp_zi) {
        // curr_zi_ord starts at 1, thus we need to subtract 1 entry
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        vi *= rand_zi.at(std::make_pair(tmp_zi, el[tmp_zi - 1]));
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
      FFInt s = nums[num_eqn - 1];

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += nums[j - 1] * b;
        t = vis[i] * t + b;
      }

      result[i] = s / t / vis[i];
    }

    // Bring result in canonical form
    ff_map poly;

    for (uint32_t i = 0; i < num_eqn; ++i) {
      poly.emplace(std::make_pair(rec_degs[i], result[i]));
    }

    nums.clear();
    return poly;
  }

  void PolyReconst::generate_anchor_points() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    rand_zi.clear();

    for (uint32_t tmp_zi = 1; tmp_zi <= n; ++tmp_zi) {
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 0), 1));
      rand_zi.emplace(std::make_pair(std::make_pair(tmp_zi, 1), get_rand()));
    }
  }

  FFInt PolyReconst::get_rand_zi(uint32_t zi, uint32_t order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return rand_zi.at(std::make_pair(zi, order));
  }

  std::vector<FFInt> PolyReconst::get_rand_zi_vec(const std::vector<uint32_t>& orders) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    std::vector<FFInt> yis {};

    for (uint32_t i = 0; i < n; ++i) {
      yis.emplace_back(rand_zi.at(std::make_pair(i + 1, orders[i])));
    }

    return yis;
  }

  bool PolyReconst::is_rand_zi_empty() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    return rand_zi.empty();
  }

  void PolyReconst::reset() {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    rand_zi = ff_pair_map();
  }

  ff_map PolyReconst::construct_tmp_canonical(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) const {
    ff_map tmp {};

    //todo hier muss irgendwo ein bug sein!
    for (auto & el : construct_canonical(ai)) { // homogenize
      std::vector<uint32_t> new_deg(n);
      new_deg[zi - 1] = el.first[0];

      for (uint32_t j = 0; j < zi - 1; j++) {
        new_deg[j] = deg_vec[j];
      }

      tmp.emplace(std::make_pair(new_deg, el.second));
    }

    return tmp;
  }

  void PolyReconst::check_for_tmp_solved_degs(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) {
    ff_map tmp = construct_tmp_canonical(deg_vec, ai);

    for (auto & el : tmp) {
      int total_deg = 0;

      for (const auto & e : el.first) total_deg += e;

      if (total_deg == deg)
        solved_degs.emplace(std::make_pair(el.first, el.second));
      else
        tmp_solved_degs.emplace(std::make_pair(el.first, el.second));
    }

    if (zi > 1) {
      std::vector<std::vector<uint32_t>>::iterator it = std::find(rec_degs.begin(), rec_degs.end(), deg_vec);
      rec_degs.erase(it);
    }
  }
}
