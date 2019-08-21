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

#pragma once

#include "FFInt.hpp"
#include "RationalNumber.hpp"

#include <vector>

namespace firefly {
  /**
   *  Parses a vector from a file with a given number of maximal entries
   *  @param line a string representing the line which should be parsed
   *  @param number_of_parameters a limiting number how many entries should be parsed
   *  @return the parsed vector
   */
  std::vector<uint32_t> parse_vector(std::string& line, int number_of_parameters = -1);
  /**
   *  Parses a vector from a file with a given number of maximal entries
   *  @param line a string representing the line which should be parsed
   *  @param number_of_parameters a limiting number how many entries should be parsed
   *  @return the parsed vector
   */
  std::vector<FFInt> parse_vector(std::string& line, std::string tmp64, int number_of_parameters = -1);
  /**
   *  Parses a rational number from a file
   *  @param line the string that should be parsed to a rational number
   *  @return the parsed rational number
   */
  RationalNumber parse_rational_number(const std::string& line);
  /**
  *  Parses a prime number counter from a file
  *  @param file_name the file name
  *  @return the parsed prime number
  */
  uint32_t parse_prime_number(const std::string& file_name);
}
