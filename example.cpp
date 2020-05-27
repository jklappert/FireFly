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

#include "firefly/DenseSolver.hpp"
#include "firefly/Reconstructor.hpp"

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
  public:
    // Constructor of the derived class
    // A default constructor is sufficient if no internal variables are required.
    // In this example a ShuntingYardParser
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

    // The evaluation of the black box
    // Return a vector of FFIntTemp objects, which are the results of the black-box evaluation
    // with values inserted for the variables.
    // In this example we compute functions which are parsed from a file with a
    // ShuntingYardParser object and the determinant of a matrix.
    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      //std::vector<FFIntTemp> result;

      // Get results from parsed expressions
      std::vector<FFIntTemp> result = par.evaluate_pre(values);

      result.emplace_back(result[0] / result[3]);

      // Build the matrix mat
      mat_ff<FFIntTemp> mat = {{result[0], result[1]}, {result[2], result[3]}};

      // Permutation vector
      std::vector<int> p {};

      // Compute LU decomposition of mat
      calc_lu_decomposition(mat, p, 2);

      // Compute determinant of mat
      result.emplace_back(calc_determinant_lu(mat, p, 2));

      return result;
    }

    // This function is called from Reconstructor when the prime field changes.
    // Update the internal variables if required.
    // In this example we precompute a few things for the ShuntingYardParser in
    // the new prime field.
    inline void prime_changed() {
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
  std::string root_dir = FIREFLY_ROOT_DIR;
  // Parse the functions from "../s_y_4_v.m" with the variables x1, y, zZ, W
  ShuntingYardParser par(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});

  // Create the user defined black box
  BlackBoxUser bb(par);

  // Initialize the Reconstructor
  Reconstructor<BlackBoxUser> reconst(4 /*n_vars*/,
                                      std::thread::hardware_concurrency() /*n_threads*/,
                                      1 /*bunch size*/,
                                      bb /*black box*//*,
                                      Reconstructor<BlackBoxUser>::CHATTY *//* verbosity mode*/);

  // Enables scan for factors
  reconst.enable_factor_scan();

  // Enables a scan for a sparse shift
  reconst.enable_shift_scan();

  // Set flag for the safe mode
  //reconst.set_safe_interpolation();

  // Write the state of all reconstruction objects after each interpolation over a prime field
  // The intermediate results are stored in ./ff_save
  //reconst.set_tags();

  // Read in all saved states from the directory ./ff_save
  //reconst.resume_from_saved_state();

  // Reconstruct the black box
  reconst.reconstruct();

  // Get results after the first interpolation. Has to be called with
  // reconst.reconstruct(1) to work
  /*for(const auto& rf : reconst.get_result_ff()) {
    std::cout << rf.to_string({"x","y","z","w"}) << "\n";
  }*/

  // Get results
  /*std::vector<RationalFunction> results = reconst.get_result();

  //std::ofstream file;
  //file.open("test.m");
  //file << "{";
  //std::string str = "";

  // Print all reconstruced functions
  for (uint32_t i = 0; i < results.size(); ++i) {
    std::cout << "Function " << i + 1 << ":\n" << results[i].to_string( {"x", "y", "z", "w"}) << "\n";

    //file << str << "\n";
    //str = results[i].to_string( {"x", "y", "z", "w"}) + "\n,";
  }

  //str.pop_back();
  //file << str << "}\n";
  //file.close();
  */

  // Rewrite result in Horner form
  /*std::string f15_horner = results[14].generate_horner({"x", "y", "z", "w"});
  std::cout << "Function 15 in Horner form:\n" << f15_horner << "\n";*/

  // Resets all statics in RatReconst to start a new reconstruction
  //RatReconst::reset();

  return 0;
}
