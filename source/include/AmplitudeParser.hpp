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
    AmplitudeParser(const std::vector<std::string>& vars);
    /**
     * 
     */
    void parse_file(const std::string& amplitude_file);
    /**
     * 
     */
    void parse_string(const std::string& amplitude, const std::list<std::string>& integral_families);
    /**
     *  Parses a rational function given as a string and validates the function on request
     *  @param fun The rational function to be parsed as a string
     *  @param vars A vector which specifies the variables of the function
     *  @param validate_fun validates the function if set to true
     */
    void parse_function(const std::string& fun, const std::vector<std::string>& vars, bool validate_fun = false);
    /**
     *  Evaluates all functions for a given parameter point and returns their result.
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate_pre(const std::vector<FFIntTemp>& values) const;
  private:
    std::unordered_map<std::string, int> vars_map {}; /**< This map holds the conversion of the used variables to an integer value to map variables to parameter points */
  };

}
