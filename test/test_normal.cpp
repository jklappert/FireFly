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

#ifdef WITH_MPI
#include "MPIWorker.hpp"
#endif

namespace firefly {
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

    virtual void prime_changed() {
      par.precompute_tokens();
    }

  private:
    // Internal variables for the black box
    // In this example a ShuntingYardParser
    ShuntingYardParser par;
    int mode = 0;
  };
}

using namespace firefly;
int main() {
#ifdef WITH_MPI
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);

  int process;
  MPI_Comm_rank(MPI_COMM_WORLD, &process);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
#endif
#ifndef WITH_MPI
  INFO_MSG("Test 1 variable");
  ShuntingYardParser p_2("../../parser_test/s_y_1_v.m", {"x"});
  BlackBoxUser b_2(p_2, 2);
  Reconstructor r_2(1, 4, b_2);
  r_2.reconstruct();
  RatReconst::reset();
  INFO_MSG("1 variable passed");

  INFO_MSG("Test normal mode");
  ShuntingYardParser p_0("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_0(p_0, 0);
  Reconstructor r_0(4, 4, b_0);
  r_0.enable_scan();
  r_0.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_0_1("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_0_1(p_0, 2);
  Reconstructor r_0_1(4, 4, b_0_1);
  r_0_1.enable_scan();
  r_0_1.reconstruct();
  INFO_MSG("Normal mode passed");
#else

  if (process == master) {
    INFO_MSG("Test 1 variable");
    ShuntingYardParser p_2("../../parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_2(p_2, 2);
    Reconstructor r_2(1, std::thread::hardware_concurrency(), b_2);
    r_2.reconstruct();
    RatReconst::reset();
    INFO_MSG("1 variable passed");

    INFO_MSG("Test normal mode");
    ShuntingYardParser p_0("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_0(p_0, 0);
    Reconstructor r_0(4, std::thread::hardware_concurrency(), b_0);
    r_0.enable_scan();
    r_0.reconstruct();
    RatReconst::reset();

    ShuntingYardParser p_0_1("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_0_1(p_0, 2);
    Reconstructor r_0_1(4, std::thread::hardware_concurrency(), b_0_1);
    r_0_1.enable_scan();
    r_0_1.reconstruct();
    INFO_MSG("Normal mode passed");
  } else {
    ShuntingYardParser p_2("../../parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_2(p_2, 2);
    MPIWorker(1, std::thread::hardware_concurrency(), b_2);
    RatReconst::reset();

    ShuntingYardParser p_0("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_0(p_0, 0);
    MPIWorker(4, std::thread::hardware_concurrency(), b_0);
    RatReconst::reset();

    ShuntingYardParser p_0_1("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_0_1(p_0, 2);
    MPIWorker(4, std::thread::hardware_concurrency(), b_0_1);
  }

  MPI_Finalize();
#endif

  return 0;
}
