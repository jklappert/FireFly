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

#include <getopt.h>
#include <string>
#include <vector>

#include "firefly/BlackBoxBase.hpp"
#include "firefly/RationalFunction.hpp"
#include "firefly/Reconstructor.hpp"

static struct option long_options[] = {
  {"parallel", required_argument, 0, 'p'},
  {"variables", required_argument, 0, 'v'},
  {0, 0, 0, 0} // mark end of array
};

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

std::vector<std::vector<std::uint64_t>> load_points(const std::string& file) {
  std::ifstream ifile(file.c_str());
  std::vector<std::vector<std::uint64_t>> points {};

  if (ifile.is_open()) {
    std::string line;

    while (std::getline(ifile, line)) {
      std::size_t pos = 0;
      int i = 0;
      std::string delimiter = " ";
      std::vector<uint64_t> tmp {};

      if (line.back() != ' ') {
        line.append(" ");
      }

      while ((pos = line.find(delimiter)) != std::string::npos) {
        tmp.emplace_back(std::stoul(line.substr(0, pos)));
        line.erase(0, pos + 1);
        ++i;
      }

      points.emplace_back(tmp);
    }
  }

  ifile.close();

  return points;
}

using namespace firefly;
int main(int argc, char* argv[]) {
  std::ofstream logger;
  logger.open("firefly-exe.log");

  std::uint32_t n = 0;
  std::uint32_t n_threads = 1;

  while (true) {
    int c;
    int option_index = 0;

    c = getopt_long_only(argc, argv, "p:v:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) {
      break;
    }

    switch (c) {
      case 'p': {
        n_threads = std::stoul(optarg);
        if (n_threads > 0) {
          INFO_MSG("Starting with " + std::to_string(n_threads) + " threads.");
          logger << "Starting with " + std::to_string(n_threads) + " threads.\n";
        } else {
          ERROR_MSG("Wrong argument for -p!");
          logger << "Wrong argument for -p!\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      case 'v': {
        n = std::stoul(optarg);
        if (n > 0) {
          INFO_MSG("Setting " + std::to_string(n) + " variables.");
          logger << "Setting " + std::to_string(n) + " variables.\n";
        } else {
          ERROR_MSG("Wrong argument for -v!");
          logger << "Wrong argument for -v!\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }
        break;
      }
    }
  }

  if (n == 0) {
    ERROR_MSG("The number of variables was not set! Use the option -v.");
    logger << "The number of variables was not set! Use the option -v.\n";
    logger.close();
    std::exit(EXIT_FAILURE);
  }

  BlackBoxFireFly bb;

  // Initialize the Reconstructor
  Reconstructor<BlackBoxFireFly> reconst(n, n_threads, 1, bb /*, Reconstructor<BlackBoxFireFly>::CHATTY */);

  // Enables scan for factors
  //reconst.enable_factor_scan();

  // Enables a scan for a sparse shift
  //reconst.enable_shift_scan();

  // Set flag for the safe mode
  //reconst.set_safe_interpolation();

  std::vector<std::vector<std::uint64_t>> anchor_points = load_points("anchor_points");
  std::vector<std::vector<std::uint64_t>> shifts = load_points("shifts");

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

  logger.close();

  return 0;
}
