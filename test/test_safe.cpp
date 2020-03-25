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

#include "DenseSolver.hpp"
#include "Reconstructor.hpp"
#include "ShuntingYardParser.hpp"
#ifdef WITH_MPI
#include "MPIWorker.hpp"
#endif

namespace firefly {
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

    virtual void prime_changed() {
      par.precompute_tokens();
    }

  private:
    ShuntingYardParser par;
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
  INFO_MSG("Test safe mode");
  ShuntingYardParser p_1("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_1(p_1);
  Reconstructor<BlackBoxUser> r_1(4, 4, b_1);
  r_1.set_safe_interpolation();
  r_1.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_1_2("../../parser_test/s_y_1_v.m", {"x"});
  BlackBoxUser b_1_2(p_1_2);
  Reconstructor<BlackBoxUser> r_1_2(1, 4, b_1_2);
  r_1_2.set_safe_interpolation();
  r_1_2.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_1_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_1_3(p_1_3);
  Reconstructor<BlackBoxUser> r_1_3(4, 4, b_1_3);
  r_1_3.set_safe_interpolation();
  r_1_3.reconstruct();
  INFO_MSG("Safe mode passed");
#else

  if (process == master) {
    INFO_MSG("Test safe mode");
    ShuntingYardParser p_1("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1(p_1);
    Reconstructor<BlackBoxUser> r_1(4, std::thread::hardware_concurrency(), b_1);
    r_1.set_safe_interpolation();
    r_1.reconstruct();
    RatReconst::reset();

    ShuntingYardParser p_1_2("../../parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_1_2(p_1_2);
    Reconstructor<BlackBoxUser> r_1_2(1, std::thread::hardware_concurrency(), b_1_2);
    r_1_2.set_safe_interpolation();
    r_1_2.reconstruct();
    RatReconst::reset();

    ShuntingYardParser p_1_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1_3(p_1_3);
    Reconstructor<BlackBoxUser> r_1_3(4, std::thread::hardware_concurrency(), b_1_3);
    r_1_3.set_safe_interpolation();
    r_1_3.reconstruct();
    INFO_MSG("Safe mode passed");
  } else {
    ShuntingYardParser p_1("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1(p_1);
    MPIWorker<BlackBoxUser>(4, std::thread::hardware_concurrency(), b_1);
    RatReconst::reset();

    ShuntingYardParser p_1_2("../../parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_1_2(p_1_2);
    MPIWorker<BlackBoxUser>(1, std::thread::hardware_concurrency(), b_1_2);
    RatReconst::reset();

    ShuntingYardParser p_1_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1_3(p_1_3);
    MPIWorker<BlackBoxUser>(4, std::thread::hardware_concurrency(), b_1_3);
  }

  MPI_Finalize();
#endif
  return 0;
}
