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

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase {
  public:
    // Constructor of the derived class
    // A default constructor is sufficient if no internal variables are required.
    // In this example a ShuntingYardParser
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

    // The evaluation of the black box
    // Return a vector of FFInt objects, which are the results of the black-box evaluation
    // with values inserted for the variables. The orderings of both vectors should
    // be fixed for all evaluations.
    // In this example we compute functions which are parsed from a file with a
    // ShuntingYardParser object and the determinant of a matrix.
    virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values) {
      //std::vector<FFInt> result;

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

    // Example for bunched evaluations
    virtual std::vector<std::vector<FFInt>> operator()(const std::vector<std::vector<FFInt>>& values) {
      // Get results from parsed expressions
      std::vector<std::vector<FFInt>> result = par.evaluate_pre(values);

      for (size_t i = 0; i != values.size(); ++i) {
        result[i].emplace_back(result[i][0] / result[i][3]);

        // Build the matrix mat
        mat_ff mat = {{result[i][0], result[i][1]}, {result[i][2], result[i][3]}};
        std::vector<int> p {};
        // Compute LU decomposition of mat
        calc_lu_decomposition(mat, p, 2);
        // Compute determinant of mat
        result[i].emplace_back(calc_determinant_lu(mat, p, 2));
      }

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
  // Parse the functions from "../s_y_test.m" with the variables x, y, z, w
  ShuntingYardParser par("../parser_test/s_y_test.m", {"x", "y", "z", "w"});
  // Parse the functions from "../s_y_1_v_test.m" with the variables x
  //ShuntingYardParser par("../parser_test/s_y_1_v_test.m", {"x"});

  // Create the user defined black box
  BlackBoxUser bb(par);

  // Initialize the Reconstructor
  Reconstructor reconst(4, 4, 1, bb/*, Reconstructor::CHATTY*/);
  // Enables a scan for a sparse shift
  reconst.enable_scan();
  // Set the safe mode
  //reconst.set_safe_interpolation();

  // Write the state of all reconstruction objects after each interpolation over a prime field
  //reconst.set_tags();

  // Give the paths to the intermediate results
  /*std::vector<std::string> file_paths = {"ff_save/0_3.txt","ff_save/1_3.txt","ff_save/2_3.txt"
    ,"ff_save/3_3.txt","ff_save/4_3.txt","ff_save/5_3.txt","ff_save/6_3.txt","ff_save/7_3.txt"};*/
  // Enables to resume from a saved state
  //reconst.resume_from_saved_state(file_paths);

  // Reconstruct the black box
  reconst.reconstruct();
  // Get results
  /*std::vector<RationalFunction> results = reconst.get_result();

  for (uint32_t i = 0; i < results.size(); ++i) {
    if(i == 4 || i == 5)
      continue;
    std::cout << "Function " << i + 1 << ":\n" << results[i].to_string( {"x", "y", "z", "w"}) << "\n";
  }

  // Rewrite result in Horner form
  std::string f8_horner = results[7].generate_horner({"x", "y", "z", "w"});
  std::cout << "Function 8 in Horner form:\n" << f8_horner << "\n";*/

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
  return 0;
}
