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

#include "AmplitudeParser.hpp"

#include <chrono>
#include <fstream>

namespace firefly {

  AmplitudeParser::AmplitudeParser() {}

  AmplitudeParser::AmplitudeParser(const std::vector<std::string>& vars) {
    for (uint32_t i = 0; i != vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
    }
  }

  void AmplitudeParser::parse_file(const std::string& amplitude_file) {
    // Check if file exists
    std::ifstream infile(amplitude_file);

    if (!infile.good()) {
      ERROR_MSG("File '" + amplitude_file + "' does not exist!");
      std::exit(EXIT_FAILURE);
    }

    std::ifstream istream;
    istream.open(amplitude_file);

    std::string line;

    while (std::getline(istream, line, ';')) {
      //parse_string(line);
    }
  };

  void AmplitudeParser::parse_string(const std::string& amplitude, const std::list<std::string>& integral_families) {
    std::string amplitude_c = amplitude;
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), ' '), amplitude_c.end());
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), '\n'), amplitude_c.end());
    std::size_t amplitude_size = amplitude_c.size();

    if (amplitude_c.length() > 0) {
      std::size_t old_pos = 0;
      std::size_t found = amplitude_c.find(']');

      while (found != std::string::npos) {
        std::string seed = amplitude_c.substr(old_pos, found + 1 - old_pos);
        std::size_t fo = seed.find('[');
        std::string int_fam = seed.substr(0, fo);

        bool found_integral = false;

        for (const auto& fam : integral_families) {
          size_t fo2 = int_fam.find(fam);
          std::string coefficient;

          if (fo2 != std::string::npos) {
            found_integral = true;
            std::string prefac = "";

            std::cout << "Integral: " << int_fam.substr(fo2) << seed.substr(fo) << "\n";

            if (fo2 != 0) {
              prefac = int_fam.substr(fo2 - 1, 1);
              std::cout << "Prefactor: " << prefac << "\n";
            }

            if (prefac == "*")
              coefficient = int_fam.substr(0, int_fam.size() - fo2 - 2);
            else if (found + 1 < amplitude_size - 1 && amplitude_c.substr(found + 1, 1) == "*") {
              std::size_t coef_end_pos = 0;

              for (const auto& tmp_fam : integral_families) {
                std::size_t tmp_found = amplitude_c.find(tmp_fam, found + 1);

                if (tmp_found != std::string::npos) {
                  coef_end_pos = tmp_found - 2;
                  break;
                }
              }

              // If this happens, we are at the end of the amplitude
              if (coef_end_pos == 0)
                coefficient = amplitude_c.substr(found + 2);
              else {
                coefficient = amplitude_c.substr(found + 2, coef_end_pos - (found + 1));
                found = coef_end_pos;
              }
            } else
              coefficient = "1";

            std::cout << "Coefficient: " << coefficient << "\n";
          }

          if (found_integral)
            break;
        }

        if (!found_integral) {
          ERROR_MSG("Unknown integral family: " + int_fam);
          std::exit(EXIT_FAILURE);
        }

        old_pos = found + 1;
        found = amplitude_c.find(']', found + 1);
      }
    } else {
      ERROR_MSG("Amplitude has no content");
      std::exit(EXIT_FAILURE);
    }
  };
}
