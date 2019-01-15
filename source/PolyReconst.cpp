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

  PolyReconst::PolyReconst(uint n_, const int deg_inp, const bool with_rat_reconst_inp) {
    type = POLY;
    n = n_;
    combined_prime = FFInt::p;
    curr_zi_order = std::vector<uint>(n, 1);

    deg = deg_inp;
    with_rat_reconst = with_rat_reconst_inp;

    for (uint i = 1; i <= n; i++) {
      std::vector<FFInt> yi;
      std::vector<PolynomialFF> ai;
      ai.reserve(300);
      ais.emplace(std::make_pair(i, std::move(ai)));
      max_deg.insert(std::make_pair(i, -1));
    }
  }

  PolyReconst::PolyReconst() {}

  void PolyReconst::set_anchor_points(const std::vector<FFInt>& anchor_points, bool force) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);

    if (rand_zi.empty() || force) {
      rand_zi.clear();

      for (uint i = 1; i <= n; i ++) {
        rand_zi.emplace(std::make_pair(std::make_pair(i, 0), 1));
        rand_zi.emplace(std::make_pair(std::make_pair(i, 1), anchor_points[i - 1]));
      }
    }
  }

  void PolyReconst::feed(const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_status);

    if (fed_prime == prime_number)
      queue.emplace_back(std::make_tuple(num, feed_zi_ord));
  }

  void PolyReconst::feed(const std::vector<FFInt>& new_yis, const FFInt& num) {
    for (uint j = 0; j < n; j++) {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      rand_zi.emplace(std::make_pair(j + 1, curr_zi_order[j]), new_yis[j]);
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

  void PolyReconst::interpolate(const FFInt& num, const std::vector<uint>& feed_zi_ord) {
    if (!done) {
      if (new_prime) {
        // be sure that you have called generate_anchor_points!
        bool runtest = true;

        solved_degs.clear();

        for (uint i = 1; i <= n; i++) {
          ais[i].clear();
        }

        for (const auto ci : combined_ci) {
          mpz_class a = ci.second;

          try {
            gi.insert(std::make_pair(ci.first, get_rational_coef(a, combined_prime)));
          } catch (const std::exception&) {
            runtest = false;
            break;
          }
        }

        if (runtest) {
          done = test_guess(num);

          if (done) {
            ais.clear();
            combined_prime = 0;
            combined_ci.clear();
            max_deg.clear();
            new_prime = false;
            use_chinese_remainder = false;
            return;
          }
        }

        gi.clear();
        zi = 1;

        if (!use_chinese_remainder) use_chinese_remainder = true;

        new_prime = false;
      }

      uint i = curr_zi_order[zi - 1] - 1;

      // Univariate Newton interpolation for the lowest stage.
      if (zi == 1) {
        if (i == 0) {
          std::vector<uint> zero_element(n);
          ff_map zero_map;
          zero_map.emplace(std::make_pair(std::move(zero_element), num));
          ais[zi].emplace_back(PolynomialFF(n, zero_map));
        } else {
          std::vector<uint> zero_element(n);
          ff_map zero_map;
          zero_map.emplace(std::make_pair(std::move(zero_element), num));
          ais[zi].emplace_back(comp_ai(zi, i, i, PolynomialFF(n, zero_map), ais[zi]));
        }

        curr_zi_order[zi - 1] ++;
      } else {
        // Build Vandermonde system
        FFInt res = num;

        for (const auto & el : solved_degs) {
          std::vector<uint> deg_vec = el.first;
          FFInt coef_num = el.second;

          for (uint tmp_zi = 1; tmp_zi < zi; tmp_zi++) {
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
          const uint order_save = curr_zi_order[zi - 1];
          curr_zi_order = std::vector<uint> (n, 1);

          for (uint i = 1; i < zi; i++) curr_zi_order[i - 1] = 0;

          curr_zi_order[zi - 1] = order_save + 1;
          ais[zi].emplace_back(comp_ai(zi, i, i, solve_transposed_vandermonde(), ais[zi]));
        } else {
          // increase all zi order of the lower stages by one
          for (uint tmp_zi = 1; tmp_zi < zi; tmp_zi++) {
            curr_zi_order[tmp_zi - 1] ++;
          }
        }
      }

      // if the lowest stage ai is zero, combine them into an ai for a higher stage
      // and check if we are done
      if (ais[zi].back().zero() || (deg != -1 && ais[zi].size() - 1 == (uint) deg)) {
        if (deg == -1 || ais[zi].back().zero())
          ais[zi].pop_back();

        if (n > 1) {
          // combine the current stage with the multivariate polynomial of the
          // previous stages and extract the reconstructed degrees to prepare
          // the gauss system
          // Remove all terms which are of total degree of the polynomial
          // to remove them from the next Vandermonde systems
          rec_degs.clear();
          PolynomialFF pol_ff = construct_canonical(zi, ais[zi]);
          PolynomialFF tmp_pol_ff = pol_ff;

          for (auto & el : tmp_pol_ff.coefs) {
            int total_deg = 0;

            for (auto & e : el.first) total_deg += e;

            if (total_deg == deg) {
              solved_degs.emplace(std::make_pair(el.first, el.second));
              pol_ff.coefs.erase(el.first);
            } else
              rec_degs.emplace_back(el.first);
          }

          // The monomials which have to be reconstructed have to
          // ordered in a monotonical way to utilize the Vandermonde
          // system solver
          std::sort(rec_degs.begin(), rec_degs.end());

          nums.reserve(rec_degs.size());

          if (rec_degs.size() == 0 && zi != n) {
            for (uint tmp_zi = zi + 1; tmp_zi <= n; tmp_zi++) {
              ais[tmp_zi].emplace_back(pol_ff);
            }

            zi = n;
          }

          if (zi != n) {
            zi ++;
            // save last interpolation of the lower stage as first a_0 of the
            // current stage
            ais[zi].emplace_back(comp_ai(zi, 0, 0, pol_ff, ais[zi]));
            // reset zi order
            curr_zi_order = std::vector<uint> (n, 1);

            for (uint i = 1; i < zi; i++) curr_zi_order[i - 1] = 0;

            curr_zi_order[zi - 1] = 2;
          } else
            check = true;
        } else if (zi == 1 && n == 1)
          check = true;

        if (check && zi == n) {
          curr_zi_order = std::vector<uint> (n, 1);
          mpz_map ci_tmp = convert_to_mpz(construct_canonical(n, ais[n]));

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

          new_prime = true;
          prime_number ++;
          check = false;
          return;
        }

        if (!with_rat_reconst) {
          auto key = std::make_pair(zi, curr_zi_order[zi - 1]);
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          set_new_rand(lock_statics, key);
        }

        return;
      }

      if (!with_rat_reconst) {
        auto key = std::make_pair(zi, curr_zi_order[zi - 1]);
        std::unique_lock<std::mutex> lock_statics(mutex_statics);
        set_new_rand(lock_statics, key);
      }
    }
  }

  Polynomial PolyReconst::get_result() {
    if (result.coefs.empty()) {
      if (done) {
        result = Polynomial(gi);
        result.sort();
        gi.clear();
      } else {
        PolynomialFF poly = construct_canonical(n, ais[n]);
        poly.coefs.insert(solved_degs.begin(), solved_degs.end());
        rn_map res;

        for (auto & el : poly.coefs) {
          res.emplace(std::make_pair(el.first, RationalNumber(el.second.n, 1)));
        }

        return Polynomial(res);
      }
    }

    return result;
  }

  PolynomialFF PolyReconst::get_result_ff() {
    PolynomialFF poly = construct_canonical(n, ais[n]);
    poly.coefs.insert(solved_degs.begin(), solved_degs.end());
    return poly;
  }

  PolynomialFF PolyReconst::comp_ai(const uint tmp_zi, int i, int ip,
                                    const PolynomialFF& num, std::vector<PolynomialFF>& ai) {
    if (ip == 0) return num;

    FFInt yi_i_p_1;
    FFInt yi_ip;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yi_i_p_1 = rand_zi[std::make_pair(tmp_zi, i + 1)];
      yi_ip = rand_zi[std::make_pair(tmp_zi, ip)];
    }

    return (comp_ai(tmp_zi, i, ip - 1, num, ai) - ai[ip - 1]) / (yi_i_p_1 - yi_ip); //yi[i + 1] - yi[ip - 1 + 1] the +1 in the first and the +1 in the second is due to the 0th element in the vector
  }

  PolynomialFF PolyReconst::construct_canonical(const uint tmp_zi, std::vector<PolynomialFF>& ai) {
    if (ai.size() == 1) return ai[0];

    return (ai[0] + iterate_canonical(tmp_zi, 1, ai));
  }

  PolynomialFF PolyReconst::iterate_canonical(const uint tmp_zi, uint i, std::vector<PolynomialFF>& ai) {
    FFInt yi;
    {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);
      yi = rand_zi[std::make_pair(tmp_zi, i)];
    }

    if (i < ai.size() - 1) {
      PolynomialFF poly = ai[i] + iterate_canonical(tmp_zi, i + 1, ai);
      return poly.mul(tmp_zi) + poly * (-yi); //yi[i - 1 + 1] +1 is due to 0th element
    }

    return ai[i] * (-yi) + ai[i].mul(tmp_zi); // yi[i - 1 + 1]
  }

  bool PolyReconst::test_guess(const FFInt& num) {
    ff_map gi_ffi = convert_to_ffint(gi);
    PolynomialFF gy(n, gi_ffi);
    std::vector<FFInt> chosen_yi(n);

    for (uint i = 1; i <= n; i++) {
      std::unique_lock<std::mutex> lock_statics(mutex_statics);

      chosen_yi[i - 1] = rand_zi[std::make_pair(i, 1)];
    }

    return gy.calc(chosen_yi) == num;
  }

  mpz_map PolyReconst::convert_to_mpz(const PolynomialFF& poly) const {
    mpz_map ci_mpz;

    for (const auto & coef : poly.coefs) {
      ci_mpz.insert(std::make_pair(coef.first, mpz_class(coef.second.n)));
    }

    return ci_mpz;
  }

  ff_map PolyReconst::convert_to_ffint(const rn_map& ri) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % FFInt::p);

      if (tmp < 0) tmp = tmp + FFInt::p;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % FFInt::p;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.insert(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  // Solves the Vandermonde linear system V*x=a
  // V is build from vis, x contain our coefficients, and a is the numerical
  // value of the function which should be interpolated for a given numerical
  // input
  PolynomialFF PolyReconst::solve_transposed_vandermonde() {
    uint num_eqn = rec_degs.size();
    std::vector<FFInt> result(num_eqn);

    if (num_eqn == 1)
      result[0] = nums[0];
    else {
      // calculate base entries of Vandermonde matrix
      std::vector<FFInt> vis;
      vis.reserve(num_eqn);

      for (const auto & el : rec_degs) {
        FFInt vi = 1;

        for (uint tmp_zi = 1; tmp_zi < zi; tmp_zi++) {
          // curr_zi_ord starts at 1, thus we need to subtract 1 entry
          std::unique_lock<std::mutex> lock_statics(mutex_statics);
          vi *= rand_zi[std::make_pair(tmp_zi, el[tmp_zi - 1])];
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

        result[i] = s / t;
      }
    }

    // Bring result in canonical form
    ff_map poly;

    for (uint i = 0; i < num_eqn; i ++) {
      poly.emplace(std::make_pair(rec_degs[i], result[i]));
    }

    nums.clear();
    return PolynomialFF(n, poly);
  }

  void PolyReconst::generate_anchor_points(uint max_order) {
    std::unique_lock<std::mutex> lock_statics(mutex_statics);
    gen_anchor_points(lock_statics, max_order);
  }
}

