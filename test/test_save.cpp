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
    BlackBoxUser(const ShuntingYardParser& par_, int mode_) : par(par_), mode(mode_) {};

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
      result.emplace_back(singular_solver(values));
      result.emplace_back(n_eq_1(values[0]));
      result.emplace_back(n_eq_4(values));
      result.emplace_back(gghh(values));
      result.emplace_back(pol_n_eq_3(values));
      result.emplace_back(ggh(values));

      return result;
    }

    virtual void prime_changed() {
      par.precompute_tokens();
      c++;

      if ((mode == 4 && c == 1) || (mode == 5 && c == 2))
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
  try {
    INFO_MSG("Test saving states and starting from them in prime 1");
    ShuntingYardParser p_4("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 4);
    Reconstructor r_4(4, 4, b_4);
    r_4.enable_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_5("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_5(p_5, 6);
    Reconstructor r_5(4, 4, b_5);
    r_5.set_tags();
    r_5.resume_from_saved_state();
    r_5.reconstruct();
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    INFO_MSG("Starting from saved states passed");
  }

  // Remove files
  remove_sates();
  remove_probes();
  RatReconst::reset();

  try {
    INFO_MSG("Test saving states and starting from them in prime 2");
    ShuntingYardParser p_4("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 5);
    Reconstructor r_4(4, 4, b_4);
    r_4.enable_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_5("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_5(p_5, 6);
    Reconstructor r_5(4, 4, b_5);
    r_5.set_tags();
    r_5.resume_from_saved_state();
    r_5.reconstruct();
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    INFO_MSG("Starting from saved states passed");
    std::cout << "\n";
  }

  // Remove files
  remove_sates();
  remove_probes();

  return 0;
}
