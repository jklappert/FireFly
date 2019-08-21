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

#include "ParserUtils.hpp"

namespace firefly {
  // TODO: Remove one function?
  std::vector<uint32_t> parse_vector(std::string& line, int number_of_parameters) {
    size_t pos = 0;
    int i = 0;
    std::string delimiter = " ";
    std::vector<uint32_t> tmp {};

    if (number_of_parameters > 0)
      tmp.reserve(number_of_parameters);

    while ((pos = line.find(delimiter)) != std::string::npos) {
      tmp.emplace_back(std::stoi(line.substr(0, pos)));
      line.erase(0, pos + 1);
      ++i;

      if (i == number_of_parameters) break;
    }

    return tmp;
  }

  std::vector<FFInt> parse_vector(std::string& line, std::string tmp64, int number_of_parameters) {
    size_t pos = 0;
    int i = 0;
    std::string delimiter = " ";
    std::vector<FFInt> tmp {};

    if (number_of_parameters > 0)
      tmp.reserve(number_of_parameters);

    while ((pos = line.find(delimiter)) != std::string::npos) {
      tmp.emplace_back(std::stoul(line.substr(0, pos)));
      line.erase(0, pos + 1);
      ++i;

      if (i == number_of_parameters) break;
    }

    return tmp;
  }

  RationalNumber parse_rational_number(const std::string& line) {
    size_t pos = line.find(" ");
    return RationalNumber(mpz_class(line.substr(0, pos)), mpz_class(line.substr(pos + 1)));
  }

  uint32_t parse_prime_number(const std::string& file_name) {
    std::string reverse_file_name = file_name;
    std::reverse(reverse_file_name.begin(), reverse_file_name.end());
    reverse_file_name.erase(0, 4);
    size_t pos = reverse_file_name.find("_");
    return std::stoi(reverse_file_name.substr(0, pos));
  }
}
