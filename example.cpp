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

#include "Reconstructor.hpp"
#include "Tests.hpp"
#include "Logger.hpp"

#include "PolyReconst.hpp"
#include "Poly.hpp"

#include "utils.hpp"

using namespace firefly;

int main() {
  // Example for the automatic interface
  // Reconstructor reconst(4, 4, Reconstructor::IMPORTANT);
  // Enables a scan for a sparse shift
  // reconst.enable_scan();
  // Give the paths to the intermediate results
  //std::vector<std::string> file_paths = {"ff_save/0_3.txt","ff_save/1_2.txt","ff_save/2_3.txt","ff_save/3_4.txt","ff_save/4_1.txt","ff_save/5_2.txt"};
  //std::vector<std::string> file_paths = {"ff_save/sing_3.txt","ff_save/n1_2.txt","ff_save/n4_3.txt","ff_save/gghh_4.txt","ff_save/pol_1.txt","ff_save/ggh_2.txt"};
  // Enables to resume from a saved state
  //reconst.resume_from_saved_state(file_paths);
  // Write the state of all reconstruction objects after each interpolation over a prime field
  //reconst.set_tags();
  // Write the state of all reconstruction objects after each interpolation over a prime field to specified tags
  //std::vector<std::string> tags = {"sing","n1","n4","gghh","pol","ggh"};
  //reconst.set_tags(tags);
  // reconst.reconstruct();
  // Get results
  /*std::vector<RationalFunction> results = reconst.get_result();
  for (auto& res : results) {
    std::cout << res << "\n";
  }*/

  // Resets all statics in RatReconst to start a new reconstruction
  // RatReconst::reset();

  // Example for the reconstruction of a rational function
  // reconstruct_rational_function();

  // Example for the reconstruction of a polynomial
  reconstruct_polynomial();
  return 0;
}

// Example of how one can use the black_box function for the automatic interface
void Reconstructor::black_box(std::vector<FFInt>& result, const std::vector<FFInt>& values) {
  result.emplace_back(singular_solver(values));
  result.emplace_back(n_eq_1(values[0]));
  result.emplace_back(n_eq_4(values));
  result.emplace_back(gghh(values));
  result.emplace_back(pol_n_eq_3(values));
  result.emplace_back(ggh(values));
}

namespace firefly {

  // Example for the reconstruction of a rational function
  void reconstruct_rational_function() {
    uint32_t n = 4;
    FFInt::set_new_prime(primes()[0]);
    RatReconst rec(n);

    // One can set a tag to start from a previously saved run after an interpolation
    // over one prime field was successful
    //rec.set_tag("rf");
    // Read in a previously saved file to resume a run from this point
    //rec.start_from_saved_file("ff_save/rf_2.txt");

    std::cout << "--------------------------------------------------------------\n";
    std::cout << "Interpolating rational function\n";
    std::cout << "--------------------------------------------------------------\n";
    // Initialize all values. t_yis are the scaled y-values and yis are the
    // unshifted ones
    std::vector<FFInt> t_yis(n);
    t_yis.reserve(n - 1);
    // t is the scaling variable
    FFInt t;

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    // One can use this option to find a sparser shift
    //rec.scan_for_sparsest_shift();

    // Feed loop
    std::vector<FFInt> shift = rec.get_zi_shift_vec();
    /*bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;

    // Generate all possible combinations of shifting variables
    auto shift_vec = generate_possible_shifts(n);

    // Run this loop until a proper shift is found
    while (!found_shift) {
      while (!rec.is_done()) {
        t = rec.get_rand();

        // Add the shift to the scaling variable
        FFInt z1 = t + shift[0];

        for (uint j = 2; j <= n; ++j) {
          t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
        }

        std::vector<FFInt> yis(n);
        yis[0] = z1;

        for (uint j = 1; j < n; ++j) {
          yis[j] = t_yis[j - 1];
        }

        //FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
        //FFInt num = n_eq_1(z1); // example for n = 1
        //FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
        //FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
        //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient

        // Feed the algorithm with the current zi_order
        ++count;
        rec.feed(t, num, rec.get_zi_order(), primes_used);
        rec.interpolate();
      }

      found_shift = rec.is_shift_working();
      counter ++;

      if (first) {
        found_shift = false;
        first = false;
        counter = 0;
      }

      rec.set_zi_shift(shift_vec[counter]);
      shift = rec.get_zi_shift_vec();
    }

    rec.set_zi_shift(shift_vec[counter - 1]);
    shift = rec.get_zi_shift_vec();
    rec.accept_shift();
    std::cout << "Total numerical runs to get sparse shift: " << count << ".\n";*/

    // In this loop the whole reconstruction of a function happens
    while (!rec.is_done()) {
      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec.get_prime()) {
        if (!rec.need_shift()) {
          rec.disable_shift();
          shift = rec.get_zi_shift_vec();
        }

        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec.get_prime());

        FFInt::set_new_prime(primes()[rec.get_prime()]);
        rec.generate_anchor_points();
        kk = 0;
      }

      // Always set the scaling variable to a random value
      t = rec.get_rand();

      // Add the shift to the scaling variable
      FFInt z1 = t + shift[0];

      for (uint j = 2; j <= n; ++j) {
        t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
      }

      std::vector<FFInt> yis(n);
      yis[0] = z1;

      for (uint j = 1; j < n; ++j) {
        yis[j] = t_yis[j - 1];
      }

      FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
      // FFInt num = n_eq_1(z1); // example for n = 1
      // FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
      // FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
      // FFInt num = bench_3(yis);
      // FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient

      // FFInt num = pol_7(yis);

      ++kk;
      ++count;

      // Feed the algorithm with the current zi_order
      rec.feed(t, num, rec.get_zi_order(), primes_used);
      rec.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    // std::cout << rec.get_result();
    std::cout << "--------------------------------------------------------------\n";
  }

  // Example for the reconstruction of a polynomial
  void reconstruct_polynomial() {
    FFInt::set_new_prime(primes()[0]);
    uint32_t n = 20;
    PolyReconst rec_poly(n);

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    std::cout << "Interpolating polynomial\n";
    std::cout << "--------------------------------------------------------------\n";
    std::vector<FFInt> yis(n);
    rec_poly.generate_anchor_points();

    // Initialize some counters
    count = 0;
    kk = 0;
    primes_used = 0;

    // Reconstruction loop
    while (!rec_poly.is_done()) {
      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec_poly.get_prime()) {
        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec_poly.get_prime());

        FFInt::set_new_prime(primes()[rec_poly.get_prime()]);
        rec_poly.generate_anchor_points();
        kk = 0;
      }

      yis = rec_poly.get_rand_zi_vec(rec_poly.get_zi_order());

      FFInt num = pol_20_20(yis);

      ++kk;
      ++count;

      rec_poly.feed(num, rec_poly.get_zi_order(), primes_used);
      rec_poly.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec_poly.get_result();
    std::cout << "--------------------------------------------------------------\n";
  }

}
