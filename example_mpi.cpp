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
//#include "Tests.hpp"

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
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
    // TODO
    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      //std::vector<FFIntTemp> result(1, FFIntTemp(42));

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

      //std::vector<FFIntTemp> b(1);
      //solve_lu(mat, p, b, 2);

      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //calc_inverse(mat, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //solve_gauss_system(mat, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //p = {};
      //mat_ff<FFIntTemp> mat2;
      //calc_inverse_lu(mat, mat2, p, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //p = {};
      //calc_determinant_lu(mat, p, 2);

      // Some functions from Test.cpp
      //result.emplace_back(singular_solver(values));
      //result.emplace_back(n_eq_1(values[0]));
      //result.emplace_back(n_eq_4(values));
      //result.emplace_back(gghh(values));
      //result.emplace_back(pol_n_eq_3(values));
      //result.emplace_back(ggh(values));

      //std::cout << "bb size " << result.size() << "\n";

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

  // Initialization of MPI processes
  // ---------------------------------------------------------------------------
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);

  int process;
  MPI_Comm_rank(MPI_COMM_WORLD, &process);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  // ---------------------------------------------------------------------------

  // Use user defined black box. Interpolate only on master and calculate probes on slaves.

  // Parse the functions from "../s_y_4_v.m" with the variables x1, y, zZ, W
  ShuntingYardParser par("../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});

  // Create the user defined black box
  BlackBoxUser bb(par);

  uint32_t bunch_size = 16;

  if (process == master) {
    Reconstructor<BlackBoxUser> reconst(4 /*n_vars*/,
                          std::thread::hardware_concurrency() /*n_threads*/,
                          bunch_size /*bunch size*/,
                          bb /*black box*//*, Reconstructor::CHATTY*/);

    reconst.enable_scan();

    reconst.reconstruct();
  } else {
    MPIWorker<BlackBoxUser>(4 /*n_vars*/,
              std::thread::hardware_concurrency() /*n_threads*/,
              bunch_size,
              bb /*black box*/);
  }

  // Finish MPI environment
  MPI_Finalize();

  return 0;
}
