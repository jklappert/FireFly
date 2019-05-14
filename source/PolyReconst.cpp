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

#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include "utils.hpp"
#include <chrono>
#include "Poly.hpp"


namespace firefly {
  // TODO check if this interpolates in combination with RatReconst to use the
  // static rand_zi of RatReconst to save additional memory -> note that
  // zi order in RatReconst and PolyReconst is not the same!
  // TODO for new prime just use Vandermonde matrices to solve interpolation problem

  ff_pair_map PolyReconst::rand_zi;
  std::mutex PolyReconst::mutex_statics;
  size_t PolyReconst::BT_threshold = 1;
  bool PolyReconst::use_Newton = true;
  bool PolyReconst::use_BT = true;

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
          Nums_for_BM.clear();
          BT_Terminator.clear();
          Lambda.clear();
          B.clear();
          L.clear();
          Delta.clear();
          BM_iteration.clear();

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

          bool finished = false;

          if(use_Newton){
	    if (i == 0)
	      ais[zero_element].emplace_back(num);
	    else
	      ais[zero_element].emplace_back(comp_ai(i, i, num, ais[zero_element]));

	    if (ais[zero_element].back() == 0) {
	      combine_res = true;
	      ais[zero_element].pop_back();
	      finished = true;
	      // std::cout << "Newton interpolation finished first in the first variable.\n";
	    } else if (deg != -1 && (uint32_t) deg == i) {
	      combine_res = true;
	      finished = true;
	    }
	    // std::cout << "Newton interpolation finished first in the first variable.\n";
	  }
	  
          if(use_BT && !finished){
            if(Nums_for_BM[zero_element].size() == 0 || i == 0){
              BT_Terminator.erase(zero_element);
              B.erase(zero_element);
              L.erase(zero_element);
              Delta.erase(zero_element);
              BM_iteration.erase(zero_element);
              Nums_for_BM.erase(zero_element);
              Lambda.erase(zero_element);
              Lambda[zero_element].emplace_back(FFInt(1));
              BT_Terminator[zero_element] = 1;
              B[zero_element].emplace_back(FFInt(0));
              L[zero_element] = 0;
              Delta[zero_element] = FFInt(1);
              BM_iteration[zero_element] = 1;
            }


            Nums_for_BM[zero_element].emplace_back(num);

            finished = Berlekamp_Massey_step(zero_element);
            if(finished){
              std::pair<std::vector<FFInt>, std::vector<size_t>> roots = rootsexponents(zero_element, get_rand_zi(zi, 1));
              if(roots.first.size() == Lambda[zero_element].size() - 1){
                std::vector<FFInt> result;
                result = solve_transposed_vandermonde2(roots.first, Nums_for_BM[zero_element]);
                //combine the result if suceeded with the former interpolated polynomial
                //check for tmp solved degrees
                check_for_tmp_solved_degs_BT(zero_element, result, roots.second);
                // std::cout << "Ben-Or and Tiwari interpolation finished first in the first variable.\n";
                combine_res = true;
                ais.erase(zero_element);
                BT_Terminator.erase(zero_element);
                B.erase(zero_element);
                L.erase(zero_element);
                Delta.erase(zero_element);
                BM_iteration.erase(zero_element);
                Nums_for_BM.erase(zero_element);
                Lambda.erase(zero_element);
              }else{
                finished = false;
                BT_Terminator[zero_element] = 1;
              }
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

            for (uint32_t tmp_zi = 1; tmp_zi <= zi; ++tmp_zi) {
              // curr_zi_ord starts at 1, thus we need to subtract 1 entry
              std::unique_lock<std::mutex> lock_statics(mutex_statics);
              coef_num *= rand_zi[std::make_pair(tmp_zi, curr_zi_order[tmp_zi - 1])].pow(deg_vec[tmp_zi - 1]);
            }

            res -= coef_num;
          }

          for (const auto & el : tmp_solved_degs) {
            std::vector<uint32_t> deg_vec = el.first;
            FFInt coef_num = el.second;

            for (uint32_t tmp_zi = 1; tmp_zi <= zi; ++tmp_zi) {
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
            uint32_t not_done_counter_BM = 0;

            for (const auto & el : solve_transposed_vandermonde()) {
              std::vector<uint32_t> key = el.first;

              bool finished = false;

              if(use_Newton){
		uint32_t tmp_deg = deg;

		for (auto ele : key) {
		  tmp_deg -= ele;
		}

		FFInt tmp_ai = comp_ai(i, i, el.second, ais[key]);
		ais[key].emplace_back(tmp_ai);

		if (tmp_ai == 0) {
		  ais[key].pop_back();
		  check_for_tmp_solved_degs(key, ais[key]);
		  ais.erase(key);
		  BT_Terminator.erase(key);
		  B.erase(key);
		  L.erase(key);
		  Delta.erase(key);
		  BM_iteration.erase(key);
		  Nums_for_BM.erase(key);
		  Lambda.erase(key);
		  finished = true;
		  // std::cout << "Newton interpolation finished first.\n";
		} else if (deg != -1 && i == tmp_deg) {
		  check_for_tmp_solved_degs(key, ais[key]);
		  ais.erase(key);
		  BT_Terminator.erase(key);
		  B.erase(key);
		  L.erase(key);
		  Delta.erase(key);
		  BM_iteration.erase(key);
		  Nums_for_BM.erase(key);
		  Lambda.erase(key);
		  finished = true;
		  // std::cout << "Newton interpolation finished first.\n";
		} else {
		  ++not_done_counter;
		}
	      }
	      
              if(use_BT && !finished){

                Nums_for_BM[key].emplace_back(el.second);

                if(Nums_for_BM[key].size() == 1){
                  BT_Terminator.erase(key);
                  B.erase(key);
                  L.erase(key);
                  Delta.erase(key);
                  BM_iteration.erase(key);
                  Lambda.erase(key);
                  Lambda[key].emplace_back(FFInt(1));
                  BT_Terminator[key] = 1;
                  B[key].emplace_back(FFInt(0));
                  L[key] = 0;
                  Delta[key] = FFInt(1);
                  BM_iteration[key] = 1;
                }

                finished = Berlekamp_Massey_step(key);

                if(finished){
                  std::pair<std::vector<FFInt>, std::vector<size_t>> roots = rootsexponents(key, get_rand_zi(zi, 1));
                  if(roots.first.size() == Lambda[key].size() - 1){
                    std::vector<FFInt> result = solve_transposed_vandermonde2(roots.first, Nums_for_BM[key]);
                    //combine the result if suceeded with the former interpolated polynomial
                    //check for tmp solved degrees
                    check_for_tmp_solved_degs_BT(key, result, roots.second);
                    // std::cout << "Ben-Or and Tiwari interpolation finished first.\n";
                    BT_Terminator.erase(key);
                    B.erase(key);
                    L.erase(key);
                    Delta.erase(key);
                    BM_iteration.erase(key);
                    Nums_for_BM.erase(key);
                    Lambda.erase(key);
                    ais.erase(key);
                  }else{
                    not_done_counter_BM++;
                    finished = false;
                    BT_Terminator[key] = 1;
                  }
                }else{
                  not_done_counter_BM++;
                }
              }
            }
            if (not_done_counter == 0 && use_Newton){
              combine_res = true;
            }

            if (not_done_counter_BM == 0 && use_BT){
              combine_res = true;
            }
	    
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

              Nums_for_BM.clear();
              BT_Terminator.clear();
              B.clear();
              L.clear();
              Lambda.clear();
              Delta.clear();
              BM_iteration.clear();

              for (const auto & el : pol_ff) {
                ais[el.first].emplace_back(el.second);
                Nums_for_BM[el.first].emplace_back(el.second);
                Lambda[el.first].emplace_back(FFInt(1));
                BT_Terminator[el.first] = 1;
                B[el.first].emplace_back(FFInt(0));
                L[el.first] = 0;
                Delta[el.first] = FFInt(1);
                BM_iteration[el.first] = 1;
                Berlekamp_Massey_step(el.first);
              }

              // reset zi order
              curr_zi_order = std::vector<uint32_t> (n, 1);
              curr_zi_order[zi - 1] = 2;
            } else {
              check = true;
              ais.clear();

              Nums_for_BM.clear();
              BT_Terminator.clear();
              Lambda.clear();
              B.clear();
              L.clear();
              Delta.clear();
              BM_iteration.clear();

              for (const auto & el : pol_ff) {
                ais[el.first].emplace_back(el.second);
                Nums_for_BM[el.first].emplace_back(el.second);
                Lambda[el.first].emplace_back(FFInt(1));
                BT_Terminator[el.first] = 1;
                B[el.first].emplace_back(FFInt(0));
                L[el.first] = 0;
                Delta[el.first] = FFInt(1);
                BM_iteration[el.first] = 1;
                Berlekamp_Massey_step(el.first);
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

              Nums_for_BM.clear();
              BT_Terminator.clear();
              Lambda.clear();
              B.clear();
              L.clear();
              Delta.clear();
              BM_iteration.clear();

              result_ff = PolynomialFF(n, tmp_pol_ff).homogenize(deg);
              result_ff.n = n + 1;
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


    result = solve_transposed_vandermonde2(vis, nums);

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

    //    std::cout << "test " << PolynomialFF(n,ff_map) << "\n";
    if (ai.size() == 1 && zi == n) {
      tmp.emplace(std::make_pair(deg_vec, ai[0]));
    } else {
      for (auto & el : construct_canonical(ai)) { // homogenize
	if(el.second != 0){
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

  void PolyReconst::check_for_tmp_solved_degs(const std::vector<uint32_t>& deg_vec, const std::vector<FFInt>& ai) {
    ff_map tmp = construct_tmp_canonical(deg_vec, ai);

    //std::cout << "Writing results after finishing Newton: "<< PolynomialFF(n, tmp);

    for (auto & el : tmp) {
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

  void PolyReconst::set_BT_threshold(size_t threshold){
    BT_threshold = threshold;
  }

  bool PolyReconst::Berlekamp_Massey_step(std::vector<uint32_t> key){
    FFInt Delta_r = 0;
    for(size_t i = 0; i < Lambda[key].size(); i++){
      Delta_r += Lambda[key].at(i) * Nums_for_BM[key].at(BM_iteration[key]-i-1);
    };
    if(Delta_r == FFInt(0)){
      B[key].insert(B[key].begin(), FFInt(0));
      while(B[key].back() == FFInt(0)){
        B[key].pop_back();
      }
      BT_Terminator[key]++;
      if(BT_Terminator[key] >= BT_threshold + 1) {
        BM_iteration[key]++;
        return true;
      }else{
        BM_iteration[key]++;
        return false;
      };
    };
    if(Delta_r != FFInt(0)){
      std::vector<FFInt> B_temp (Lambda[key]);
      std::vector<FFInt> Lambda_temp;
      B[key].insert(B[key].begin(), FFInt(0));
      for(size_t j = 0; j < Lambda[key].size() || j < B[key].size(); j++){
        if(j < Lambda[key].size() && j < B[key].size()){
          Lambda_temp.emplace_back(Lambda[key].at(j) - Delta_r/Delta[key]*B[key].at(j));
        };
        if(j < Lambda[key].size() && j >= B[key].size()){
          Lambda_temp.emplace_back(Lambda[key].at(j));
        };
        if(j >= Lambda[key].size() && j < B[key].size()){
          Lambda_temp.emplace_back(- Delta_r/Delta[key]*B[key].at(j));
        };
      };
      if(2*L[key] < BM_iteration[key]){
        L[key] = BM_iteration[key] - L[key];
        Delta[key] = Delta_r;
        B[key].swap(B_temp);
      }
      Lambda[key].swap(Lambda_temp);
      while(Lambda[key].back() == FFInt(0)){
        Lambda[key].pop_back();
      };
      while(B[key].back() == FFInt(0)){
        B[key].pop_back();
      };
      BM_iteration[key]++;
      BT_Terminator[key] = 1;
      return false;
    }
  }

  std::pair<std::vector<FFInt>, std::vector<size_t>> PolyReconst::rootsexponents(std::vector<uint32_t> key, FFInt base){
    std::pair<std::vector<FFInt>, std::vector<size_t>> roots;
    FFInt a(1);
    size_t count = 0;
    size_t sol_deg = 0;
    for(size_t i = 0; i < key.size(); i++){
      sol_deg += key.at(i);
    }
    while(roots.first.size() < Lambda[key].size() - 1){
      FFInt result(0);
      for(size_t i = 0; i < Lambda[key].size(); i++){
        result += Lambda[key].at(Lambda[key].size() - i- 1) * a.pow(i);
      };
      if(result == 0){
        roots.first.emplace_back(a);
        roots.second.emplace_back(count);
        rand_zi.emplace(std::make_pair(std::make_pair(zi, count), a));
      }
      a *= base;
      if(a == FFInt(1)){break;};
      count++;
      if(count > deg- sol_deg){break;};
    }
    if(roots.first.size() != Lambda[key].size()-1){
      std::cout << "\033[1;31mThe Polynomial calculated by the Berlekamp/Massey algorithm is not correct\033[0m\n";
    }
    return roots;
  }

  std::pair<std::vector<FFInt>, std::vector<size_t>> PolyReconst::rootsexponents2(std::vector<uint32_t> key, FFInt base){
    std::pair<std::vector<FFInt>, std::vector<size_t>> rootsexponent;
    Poly LambdaPoly(Lambda[key]);
    LambdaPoly.rev();
    std::vector<FFInt> roots = LambdaPoly.roots();
    FFInt a(1);
    for(size_t count = 0; count <= deg; count++){
      for(size_t j = 0; j < roots.size(); j++){
        if(a == roots.at(j)){
          rootsexponent.first.emplace_back(a);
          rootsexponent.second.emplace_back(count);
        }
      }
      a *= base;
      if(a == FFInt(1)){break;};
      if(rootsexponent.first.size() == roots.size()){break;};
    }
    if(rootsexponent.first.size() != LambdaPoly.get_deg()){
      std::cout << "\033[1;31mThe Polynomial calculated by the Berlekamp/Massey algorithm is not correct\033[0m\n";
    }
    return rootsexponent;
  }

  std::vector<FFInt> PolyReconst::solve_transposed_vandermonde2(std::vector<FFInt> vis, std::vector<FFInt> fis){
    uint32_t num_eqn = vis.size();

    if(num_eqn == 0){
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

  void PolyReconst::check_for_tmp_solved_degs_BT(const std::vector<uint32_t> & deg_vec, const std::vector<FFInt> & coeffs, std::vector<size_t> & exponents){
    ff_map tmp;
    
    for(size_t i = 0; i < exponents.size(); i++){
      std::vector<uint32_t> new_deg_vec = deg_vec;
      new_deg_vec[zi-1] = exponents.at(i);
      tmp.emplace(std::make_pair(new_deg_vec, coeffs.at(i)));
    }
    if(exponents.size() == 0){
      tmp.emplace(std::make_pair(deg_vec, FFInt(0)));
    }
    //std::cout << "Writing results after finishing bt: "<< PolynomialFF(n, tmp);
    for (auto & el : tmp) {
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

  void PolyReconst::set_Newton(bool use_Newton_new){
    use_Newton = use_Newton_new;
  }

  void PolyReconst::set_BT(bool use_BT_new){
    use_BT = use_BT_new;
  }
}

