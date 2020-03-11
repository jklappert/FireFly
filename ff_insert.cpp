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
#include <sys/stat.h>

using namespace firefly;
int main(int argc, char *argv[]) {
  auto time0 = std::chrono::high_resolution_clock::now();
  uint32_t n_threads = 1;
  uint32_t bs = 1;
  bool factor_scan = true;
  bool save_mode = false;
  bool no_interpolation = false;

  for (int i = 0; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-p" || arg == "--parallel") {
      n_threads = std::stoi(argv[i + 1]);
    }

    if (arg == "-bs" || arg == "--bunchsize")
      bs = std::stoi(argv[i + 1]);

    if (arg == "-nfs" || arg == "--nofactorscan")
      factor_scan = false;

    if (arg == "-ni" || arg == "--nointerpolation")
      no_interpolation = true;

    if (arg == "-s" || arg == "--save")
      save_mode = true;

    if (arg == "-h" || arg == "--help") {
      std::cerr << "Usage: " << argv[0] << "\n"
                << "Options:\n  -p,--parallel Sets the number of used threads\n"
                << "  -bs,--bunchsize         Sets the maximum bunch size\n"
                << "  -nfs,--nofactorscan     Disables the factor scan\n"
                << "  -ni,--nointerpolation   Disables the interpolation and writes coefficients to files\n"
                << "  -s,--save               Enables the storage of intermediate results\n";
      return 1;
    }
  }

  // Parse integral families
  std::ifstream fam_file_test("config/functions");
  std::ifstream fam_file;


  if (!fam_file_test.good()) {
    ERROR_MSG("Please add a file definining the occurring functions in 'config/functions'");
    std::exit(EXIT_FAILURE);
  }

  std::vector<std::string> families{};

  fam_file.open("config/functions");

  std::string line;

  while (std::getline(fam_file, line, '\n')) {
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    if (line.size() != 0)
      families.emplace_back(line);
  }

  fam_file.close();

  // Parse variables
  std::ifstream var_file_test("config/vars");
  std::ifstream var_file;

  if (!var_file_test.good()) {
    ERROR_MSG("Please add a file definining the occurring variables 'config/vars'");
    std::exit(EXIT_FAILURE);
  }

  std::vector<std::string> vars{};

  var_file.open("config/vars");

  while (std::getline(var_file, line, '\n')) {
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    if (line.size() != 0)
      vars.emplace_back(line);
  }

  var_file.close();

  // Construct the amplitude parser
  AmplitudeParser ap(vars, families);

  // Parse the amplitude
  ap.parse_amplitude_file("config/in.m");

  // Parse IBP tables
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "replacements");

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
    ap.parse_ibp_table_file("replacements/" + file);
  }

  // Build the black boxes and start reconstruction
  size_t masters = ap.check_for_unreplaced_masters();

  if (!no_interpolation) {
    std::ofstream file;
    file.open("out.m");
    file << "{\n";
    file.close();

    for (size_t i = 0; i != masters; ++i) {
      if (i == 0)
        INFO_MSG("Reconstructing coefficient of basis function: " + ap.get_master(i) + "\n");

      auto bb = ap.build_black_box(i);

      // Construct the reconstructor
      Reconstructor<FFAmplitudeBlackBox> reconst(bb.n, n_threads, bs, bb);

      // Settings
      if (factor_scan)
        reconst.enable_factor_scan();

      reconst.enable_scan();

      bool renamed_ff_save = false;

      if (save_mode) {
        std::string tmp = "ff_save_" + ap.get_master(i);
        struct stat buffer;

        if (stat(tmp.c_str(), &buffer) == 0) {
          struct stat buffer_2;

          if (stat("ff_save", &buffer_2) == 0) {
            std::rename("ff_save", "ff_save_tmp");
            renamed_ff_save = true;
          }

          std::rename(tmp.c_str(), "ff_save");
        }

        reconst.set_tags({ap.get_master(i)});
        reconst.resume_from_saved_state();
      }

      // Reconstruct
      reconst.reconstruct();

      std::vector<RationalFunction> results = reconst.get_result();
      file.open("out.m", std::ios_base::app);

      if (!results.back().zero())
        file <<  "+ " << ap.get_master(i) << "*" + results.back().generate_horner(vars) << "\n";

      file.close();

      if (save_mode) {
        std::string tmp = "ff_save_" + ap.get_master(i);
        std::rename("ff_save", tmp.c_str());

        if (renamed_ff_save)
          std::rename("ff_save_tmp", "ff_save");
      }

      RatReconst::reset();

      if (i + 1 != masters) {
        std::cout << "\n";
        INFO_MSG("Reconstructing coefficient of basis function: " + ap.get_master(i + 1) + "\n");
      }

    }

    file.open("out.m", std::ios_base::app);
    file << "}\n";
    file.close();

    std::cout << "\n";
    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Reconstructed expression in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s");

    INFO_MSG("Result has been written to 'out.m'");
  } else {
    mkdir("coefficients", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (size_t i = 0; i != masters; ++i) {
      std::ofstream coef_file;
      std::string file_name = "coefficients/" + ap.get_master(i) + ".m";
      coef_file.open(file_name.c_str());
      coef_file << "{\n" << ap.get_master(i) + "*(" << ap.get_unsimplified_coef(i) << ")\n}\n";
      coef_file.close();
    }

    INFO_MSG("Unsimplified coefficients written to 'coefficients' directory");
  }

  return 0;
}
