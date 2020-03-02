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

#pragma once

#include "FFInt.hpp"
#include "Logger.hpp"

#include <unordered_map>
#include <unordered_set>
#include <list>

namespace firefly {
    /**
   * @class AmplitudeParser
   * @brief A parser for amplitudes that generates replacement rules for integrals and their coefficients which can be then interpolated over finite fields
   */
  class AmplitudeParser {
  public:
    /**
     *  Default constructor
     */
    AmplitudeParser();
    /**
     *  Constructor which parses a list of rational functions and prepares them for evaluation.
     *  @param file The path to the file in which the rational functions are stored
     *  @param vars A vector which specifies the variables of the functions. Each variable is allowed to be built of lower and upper case letters and numbers up to 16 characters.
     */
    AmplitudeParser(const std::vector<std::string>& vars, const std::list<std::string>& integral_families_);
    /**
     * 
     */
    void parse_file(const std::string& amplitude_file);
    /**
     * 
     */
    void parse_string(const std::string& amplitude);
    /**
     * 
     */
    void parse_ibp_table_file(const std::string& ibp_table);
    /**
     * 
     */
    void parse_ibp_table_string(const std::string& ibp_table);
  private:
    std::unordered_map<std::string, int> vars_map {}; /**< This map holds the conversion of the used variables to an integer value to map variables to parameter points */
    std::list<std::string> integral_families{};
  };
}
