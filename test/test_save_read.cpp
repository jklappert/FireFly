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
#include "firefly/ShuntingYardParser.hpp"
#include "firefly/tinydir.h"

#ifdef WITH_MPI
#include "firefly/MPIWorker.hpp"
#endif

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
  public:
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      //std::vector<FFInt> result;

      // Get results from parsed expressions
      std::vector<FFIntTemp> result = par.evaluate_pre(values);

      result.emplace_back(result[0] / result[3]);

      // Build the matrix mat
      mat_ff<FFIntTemp> mat = {{result[0], result[1]}, {result[2], result[3]}};
      std::vector<int> p {};
      // Compute LU decomposition of mat
      calc_lu_decomposition(mat, p, 2);
      // Compute determinant of mat
      result.emplace_back(calc_determinant_lu(mat, p, 2));

      return result;
    }

    inline void prime_changed() {
      par.precompute_tokens();
      c++;
    }

  private:
    // Internal variables for the black box
    // In this example a ShuntingYardParser
    ShuntingYardParser par;
    int c = 0;
  };
}

void remove_states() {
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

void remove_tmp() {
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "ff_save/tmp");

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
    paths.emplace_back("ff_save/tmp/" + file);
  }

  for (const auto & el : paths) {
    std::remove(el.c_str());
  }
}

using namespace firefly;
int main() {
  std::string root_dir = FIREFLY_ROOT_DIR;
#ifdef WITH_MPI
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);

  int process;
  MPI_Comm_rank(MPI_COMM_WORLD, &process);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if (process == master) {
    INFO_MSG("Test saving states and starting from them in prime 1");
    ShuntingYardParser p_4(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4);
    Reconstructor<BlackBoxUser> r_4(4, 4, b_4);
    r_4.enable_shift_scan();
    r_4.set_tags();
    r_4.resume_from_saved_state();
    r_4.reconstruct();

    // Remove files
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    remove_states();
    remove_probes();
    remove_tmp();
  } else {
    ShuntingYardParser p_1(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1(p_1);
    MPIWorker<BlackBoxUser>(4, std::thread::hardware_concurrency(), b_1);
  }

  MPI_Finalize();
#endif

  // Clean up
  std::remove("ff_save");
  std::remove("firefly.log");

  return 0;
}
