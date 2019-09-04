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
#include "tinydir.h"

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase {
  public:
    // Constructor of the derived class
    // A default constructor is sufficient if no internal variables are required.
    // In this example a ShuntingYardParser
    BlackBoxUser(const ShuntingYardParser& par_, int mode_) : par(par_), mode(mode_) {};

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
      if (mode == 0) {
        result.emplace_back(singular_solver(values));
        result.emplace_back(n_eq_1(values[0]));
        result.emplace_back(n_eq_4(values));
        result.emplace_back(gghh(values));
        result.emplace_back(pol_n_eq_3(values));
        result.emplace_back(ggh(values));
      }

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
      c++;

      if (mode == 4 && c == 1)
        throw std::runtime_error("Abort for save test.");
    }

  private:
    // Internal variables for the black box
    // In this example a ShuntingYardParser
    ShuntingYardParser par;
    int mode = 0;
    int c = 0;
  };
}

void remove_sates() {
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "ff_save/states");

  std::vector<std::string> files;
  std::vector<std::string> paths;

  for (size_t i = 0; i != dir.n_files; ++i) {
    tinydir_file file;
    tinydir_readfile_n(&dir, &file, i);

    if (!file.is_dir) {
      files.emplace_back(file.name);
    }
  }

  tinydir_close(&dir);

  for (const auto & file : files) {
    paths.emplace_back("ff_save/states/" + file);
  }

  for (const auto & el : paths) {
    std::remove(el.c_str());
  }
}

void remove_probes() {
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "ff_save/probes");

  std::vector<std::string> files;
  std::vector<std::string> paths;

  for (size_t i = 0; i != dir.n_files; ++i) {
    tinydir_file file;
    tinydir_readfile_n(&dir, &file, i);

    if (!file.is_dir) {
      files.emplace_back(file.name);
    }
  }

  tinydir_close(&dir);

  for (const auto & file : files) {
    paths.emplace_back("ff_save/probes/" + file);
  }

  for (const auto & el : paths) {
    std::remove(el.c_str());
  }
}

using namespace firefly;
int main() {
  INFO_MSG("Test safe mode");
  ShuntingYardParser p_1("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_1(p_1, 1);
  Reconstructor r_1(4, 4, b_1, Reconstructor::SILENT);
  r_1.set_safe_interpolation();
  r_1.reconstruct();
  RatReconst::reset();
  INFO_MSG("Safe mode passed");

  INFO_MSG("Test 1 variable");
  ShuntingYardParser p_2("../../parser_test/s_y_1_v.m", {"x"});
  BlackBoxUser b_2(p_2, 2);
  Reconstructor r_2(1, 4, b_2, Reconstructor::SILENT);
  r_2.reconstruct();
  RatReconst::reset();
  INFO_MSG("1 variable passed");

  INFO_MSG("Test normal mode");
  ShuntingYardParser p_0("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_0(p_0, 0);
  Reconstructor r_0(4, 4, b_0, Reconstructor::SILENT);
  r_0.enable_scan();
  r_0.reconstruct();
  RatReconst::reset();
  INFO_MSG("Normal mode passed");

  INFO_MSG("Test bunched evaluation");
  ShuntingYardParser p_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_3(p_3, 3);
  Reconstructor r_3(4, 4, 4, b_3, Reconstructor::SILENT);
  r_3.enable_scan();
  r_3.reconstruct();
  RatReconst::reset();
  ShuntingYardParser p_3_2("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_3_2(p_3_2, 3);
  Reconstructor r_3_2(4, 4, 4, b_3_2, Reconstructor::SILENT);
  r_3_2.set_safe_interpolation();
  r_3_2.reconstruct();
  RatReconst::reset();
  INFO_MSG("Bunched evaluation passed");

  try {
    INFO_MSG("Test saving states and starting from them");
    ShuntingYardParser p_4("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 4);
    Reconstructor r_4(4, 4, b_4, Reconstructor::SILENT);
    r_4.enable_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_5("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_5(p_5, 5);
    Reconstructor r_5(4, 4, b_5, Reconstructor::SILENT);
    r_5.set_tags();
    r_5.resume_from_saved_state();
    r_5.reconstruct();
    std::remove("ff_save/validation");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    INFO_MSG("Starting from saved states passed");
    std::cout << "\n";

    // Remove files
    remove_sates();

    remove_probes();

    INFO_MSG("All tests passed");
  }

  return 0;
}
