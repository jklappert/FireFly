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

#include "firefly/Reconstructor.hpp"

namespace firefly {
  // Empty black box
  class BlackBoxFireFly : public BlackBoxBase<BlackBoxFireFly> {
  public:
    BlackBoxFireFly() {};

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      std::vector<FFIntTemp> result;
      return result;
    }

    inline void prime_changed() {};
  };
}

using namespace firefly;
int main() {
  std::string root_dir = FIREFLY_ROOT_DIR;

  BlackBoxFireFly bb;

  // Initialize the Reconstructor
  Reconstructor<BlackBoxFireFly> reconst(4 /*n_vars*/,
                                      std::thread::hardware_concurrency() /*n_threads*/,
                                      1 /*bunch size*/,
                                      bb /*black box*//*,
                                      Reconstructor<BlackBoxFireFly>::CHATTY *//* verbosity mode*/);

  // Enables scan for factors
  //reconst.enable_factor_scan();

  // Enables a scan for a sparse shift
  //reconst.enable_shift_scan();

  // Set flag for the safe mode
  //reconst.set_safe_interpolation();

  std::vector<std::vector<std::uint64_t>> anchor_points {{224234, 23478923478, 2394789234}, {224234, 23478923478, 2394789234}, {224234, 23478923478, 2394789234}};
  std::vector<std::vector<std::uint64_t>> shifts {{25, 224234, 23478923478, 2394789234}, {25, 224234, 23478923478, 2394789234}, {25, 224234, 23478923478, 2394789234}};

  reconst.set_anchor_points(anchor_points);
  reconst.set_shifts(shifts);

  reconst.load_precomputed_probes();

  // Write the state of all reconstruction objects after each interpolation over a prime field
  // The intermediate results are stored in ./ff_save
  reconst.set_tags();

  // Read in all saved states from the directory ./ff_save
  reconst.resume_from_saved_state();

  reconst.load_precomputed_probes();

  // Reconstruct the black box
  reconst.reconstruct();

  // Print all found factors
  /*for (const auto & el : reconst.get_factors_string({"x","y","z","w"})) {
    std::cout << el << "\n";
  }*/

  // Get results after the first interpolation. Has to be called with
  // reconst.reconstruct(1) to work
  /*for(const auto& rf : reconst.get_result_ff()) {
    std::cout << rf.to_string({"x","y","z","w"}) << "\n";
  }*/

  // Get results
  std::vector<RationalFunction> results = reconst.get_result();

  //std::ofstream file;
  //file.open("test.m");
  //file << "{";
  //std::string str = "";

  // Print all reconstruced functions
  for (uint32_t i = 0; i != results.size(); ++i) {
    std::cout << "Function " << i + 1 << ":\n" << results[i].to_string( {"x", "y", "z", "w"}) << "\n";

    //file << str << "\n";
    //str = results[i].to_string( {"x", "y", "z", "w"}) + "\n,";
  }

  //str.pop_back();
  //file << str << "}\n";
  //file.close();

  // Rewrite result in Horner form
  /*std::string f15_horner = results[14].generate_horner({"x", "y", "z", "w"});
  std::cout << "Function 15 in Horner form:\n" << f15_horner << "\n";*/

  // Resets all statics in RatReconst to start a new reconstruction
  //RatReconst::reset();

  return 0;
}
