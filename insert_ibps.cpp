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

#include "Reconstructor.hpp"
#include "AmplitudeParser.hpp"

using namespace firefly;
int main(int argc, char *argv[]) {
  uint32_t n_threads = 1;
  uint32_t bs = 1;
  bool factor_scan = true;

  for (int i = 0; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-p" || arg == "--parallel") {
      n_threads = std::stoi(argv[i + 1]);
    }

    if (arg == "-bs" || arg == "--bunchsize")
      bs = std::stoi(argv[i + 1]);

    if (arg == "-nfs" || arg == "--nofactorscan")
      factor_scan = false;

    if (arg == "-h" || arg == "--help") {
      std::cerr << "Usage: " << argv[0] << "\n"
                << "Options:\n  -p,--parallel Sets the number of used threads\n"
                << "  -bs,--bunchsize         Sets the maximum bunch size\n"
                << "  -nfs,--nofactorscan     Disables the factor scan\n";
      return 1;
    }
  }

  // Parse integral families
  std::ifstream fam_file_test("amplitude/families");
  std::ifstream fam_file;


  if (!fam_file_test.good()) {
    ERROR_MSG("Please add a file definining the occurring integral families in 'amplitude/families'");
    std::exit(EXIT_FAILURE);
  }

  std::vector<std::string> families{};

  fam_file.open("amplitude/families");

  std::string line;

  while (std::getline(fam_file, line, '\n')) {
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    if (line.size() != 0)
      families.emplace_back(line);
  }

  fam_file.close();

  // Parse variables
  std::ifstream var_file_test("amplitude/vars");
  std::ifstream var_file;

  if (!var_file_test.good()) {
    ERROR_MSG("Please add a file definining the occurring variables 'amplitude/vars'");
    std::exit(EXIT_FAILURE);
  }

  std::vector<std::string> vars{};

  var_file.open("amplitude/vars");

  while (std::getline(var_file, line, '\n')) {
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    if (line.size() != 0)
      vars.emplace_back(line);
  }

  var_file.close();

  // Construct the amplitude parser
  AmplitudeParser ap(vars, families);

  // Parse the amplitude
  ap.parse_amplitude_file("amplitude/amplitude_in.m");

  // Parse IBP tables
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "amplitude/ibps");

  std::vector<std::string> files;

  for (size_t i = 0; i != dir.n_files; ++i) {
    tinydir_file file;
    tinydir_readfile_n(&dir, &file, i);

    if (!file.is_dir) {
      files.emplace_back(file.name);
    }
  }

  tinydir_close(&dir);

  for (const auto & file : files) {
    ap.parse_ibp_table_file("amplitude/ibps/" + file);
  }

  // Build the black box
  auto bb = ap.build_black_box();

  // Construct the reconstructor
  Reconstructor<FFAmplitudeBlackBox> reconst(bb.n, n_threads, bs, bb);

  // Settings
  if (factor_scan)
    reconst.enable_factor_scan();

  reconst.enable_scan();

  // Reconstruct
  reconst.reconstruct();

  std::ofstream file;
  file.open("amplitude_out.m");
  std::vector<RationalFunction> results = reconst.get_result();
  file << "{\n";

  if (results.size() == 0) {
    file << "0\n";
  } else {
    for (uint32_t i = 0; i < results.size(); ++i) {
      file <<  "+ " << ap.get_master(i) << "*" + results[i].to_string(vars) << "\n";
    }
  }

  file << "}\n";
  file.close();
  INFO_MSG("Result has been written to 'amplitude_out.m'");

  return 0;
}
