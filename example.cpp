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

#include "DenseSolver.hpp"
#include "Reconstructor.hpp"
#include "ShuntingYardParser.hpp"
#include "Tests.hpp"
#include "utils.hpp"

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase {
  public:
    // Constructor of the derived class
    // A default constructor is sufficient if no internal variables are required.
    // In this example a ShuntingYardParser
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

    // The evaluation of the black box
    // Return a vector of FFInt, which are the results of the black-box evaluation
    // with values inserted for the variables. The orderings of both vectors should
    // be fixed for all evaluations.
    // In this example we compute functions which are parsed from a file with a
    // ShuntingYardParser and the determinant of a matrix.
    virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values) {
      //std::vector<FFInt> result {};

      // Get results from parsed expressions
      std::vector<FFInt> result = par.evaluate_pre(values);
      result.emplace_back(result[0] / result[3]);

      // Build the matrix mat
      mat_ff mat = {{result[0], result[1]}, {result[2], result[3]}};
      std::vector<int> p {};
      // Compute LU decomposition of mat
      calc_lu_decomposition(mat, p, 2);
      // Compute determinant of mat
      result.emplace_back(calc_determinant_lu(mat, p, 2));

      // Some functions from Test.cpp
      /*result.emplace_back(singular_solver(values));
      result.emplace_back(n_eq_1(values[0]));
      result.emplace_back(n_eq_4(values));
      result.emplace_back(gghh(values));
      result.emplace_back(pol_n_eq_3(values));
      result.emplace_back(ggh(values));*/

      return result;
    }

    // This function is called from Reconstructor when the prime field changes.
    // Update the internal variables if required.
    // In this example we precompute a few things for the ShuntingYardParser in
    // the new prime field.
    virtual void prime_changed() {
      par.precompute_tokens();
    }

  private:
    // Internal variables for the black box
    // In this example a ShuntingYardParser
    ShuntingYardParser par;
  };
}

using namespace firefly;
int main() {
  // Example of ShuntingYardParser
  // Parse the functions from "../s_y_test.m" with the variables x, y, z
  ShuntingYardParser par("../s_y_test.m", {"x", "y", "z"});

  // Create the user defined black box
  //BlackBoxUser bb(par);

  // Initialize the Reconstructor
  //Reconstructor reconst(3, 4, bb/*, Reconstructor::CHATTY*/);
  // Enables a scan for a sparse shift
  //reconst.enable_scan();
  // Set the safe mode
  //reconst.set_safe_interpolation();

  // Give the paths to the intermediate results
  //std::vector<std::string> file_paths = {"ff_save/0_0.txt","ff_save/1_0.txt","ff_save/2_0.txt","ff_save/3_0.txt","ff_save/4_0.txt","ff_save/5_0.txt"};
  // Enables to resume from a saved state
  //reconst.resume_from_saved_state(file_paths);
  // Write the state of all reconstruction objects after each interpolation over a prime field
  //reconst.set_tags();

  // Reconstruct the black box
  //reconst.reconstruct();
  // Get results
  /*std::vector<RationalFunction> results = reconst.get_result();

  for (uint32_t i = 0; i < results.size(); ++i) {
    std::cout << "Function " << i + 1 << ":\n" << results[i].to_string( {"x", "y", "z"}) << "\n";
  }*/

  // Rewrite result in Horner form
  //std::string f6_horner = results[5].generate_horner({"x","y","z"});
  //std::cout << "Function 6 in Horner form:\n" << f6_horner << "\n";

  // Resets all statics in RatReconst to start a new reconstruction
  //RatReconst::reset();

  // Some examples of the dense solver functions

  // Initialize matrices
  /*mat_ff a = {{0*33,17,25},{59595, 989983749,99},{23213213, 4354354353,0*434232324}};
  mat_ff a_2 = a;
  mat_ff a_3 = a;
  a = {{0*33,17,25,10},{59595, 989983749,99,14},{23213213, 4354354353,0*434232324,200}};
  mat_ff inv;
  //----------------------------------------------------------------------------
  // Calculate inverse of a using Gauss-Jordan algorithm
  calc_inverse(a_3, 3);
  // Print inverse of a_3 obtained with Gauss-Jordan
  std::cout << "Inverse\n" << a_3[0][0] << " " << a_3[0][1] << " " << a_3[0][2] << "\n"
  << a_3[1][0] << " " << a_3[1][1] << " " << a_3[1][2] << "\n"
  << a_3[2][0] << " " << a_3[2][2] << " " << a_3[2][2] << "\n";
  // Solve a using Gauss algorithm
  std::vector<FFInt> res = solve_gauss_system(a, 3);
  std::cout << "Sol Gauss\n" << res[0] << " " << res[1] << " " << res[2] << "\n";
  //----------------------------------------------------------------------------
  // Initialize permutation matrix for LU decompositions
  std::vector<int> p;
  // Initialize vector for solution of LU system
  std::vector<FFInt> b = {10,14,200};
  // Calculate LU decomposition of a_2
  calc_lu_decomposition(a_2, p, 3);
  // Solve a_2 using LU decomposition
  res = solve_lu(a_2, p, b, 3);
  std::cout << "Sol LU\n" << res[0] << " " << res[1] << " " << res[2] << "\n";
  // Calculate inverse of a_2 using LU decomposition
  calc_inverse_lu(a_2, inv, p, 3);
  // Calculate determinant of a_2 using LU decomposition
  std::cout << "Det LU " << calc_determinant_lu(a_2, p, 3) << "\n";
  // Print inverse obtained using LU decomposition
  std::cout << "Inverse LU\n" << inv[0][0] << " " << inv[0][1] << " " << inv[0][2] << "\n"
  << inv[1][0] << " " << inv[1][1] << " " << inv[1][2] << "\n"
  << inv[2][0] << " " << inv[2][1] << " " << inv[2][2] << "\n";*/
  reconstruct_rational_function();
  return 0;
}

namespace firefly{
  void reconstruct_rational_function() {
    uint32_t n = 3;
    FFInt::set_new_prime(primes()[0]);
    BaseReconst br;
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    br.set_seed(seed);
    ShuntingYardParser par("../s_y_test.m", {"x", "y", "z"});

    RatReconst rec(n);
    //rec.set_safe_interpolation();

    // One can set a tag to start from a previously saved run after an interpolation
    // over one prime field was successful
    //rec.set_tag("");
    // Read in a previously saved file to resume a run from this point
    //rec.start_from_saved_file("ff_save/1_1.txt");

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
        FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
        //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient
        //FFInt num = (FFInt(primes()[0]) + FFInt(primes()[2]) * FFInt(primes()[1]) * yis[0] + 3 * yis[0].pow(2) + yis[1]) / (yis[1]);

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

        par.precompute_tokens();
      }

      // Always set the scaling variable to a random value
      t = rec.get_rand();

      // Add the shift to the scaling variable
      FFInt z1 = t + shift[5];

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
      //FFInt num = bench_3(yis);
      //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient
      //FFInt num = FFInt(primes()[1]) * FFInt(primes()[3]) * (FFInt(primes()[0]) + FFInt(primes()[2]) * FFInt(primes()[1]) * yis[0] + 3 * yis[0].pow(2) + yis[1]) / (FFInt(primes()[1]) + yis[1]);
      FFInt num = par.evaluate_pre(yis)[1];
      ++kk;
      ++count;

      // Feed the algorithm with the current zi_order
      rec.feed(t, num, rec.get_zi_order(), primes_used);
      rec.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec.get_result();
    std::cout << "--------------------------------------------------------------\n";
  }
}
