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

#include "firefly/BlackBoxBase.hpp"
#include "firefly/FFInt.hpp"
#include "firefly/Logger.hpp"
#include "firefly/ShuntingYardParser.hpp"

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <list>

namespace firefly {
  namespace coef_type {
    enum coef_type {
      PREFACTOR = 0,
      REPEATED_REP = 1,
      COEF_TYPE_SIZE = 2,
    }; /**< The key to identify a prefactor from a MI coefficient */
  }

  /**
   * @class FFAmplitudeBlackBox
   * @brief The black box of the amplitude parser which inserts IBP relations into an amplitude and returns the reduced amplitude in terms of master integrals
   */
  class FFAmplitudeBlackBox : public BlackBoxBase<FFAmplitudeBlackBox> {
  public:
    FFAmplitudeBlackBox(const std::vector<std::string>& vars,
                        const std::vector<std::string>& functions) {
      par = ShuntingYardParser(functions, vars);
      n = vars.size();
    };

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      return par.evaluate_pre(values);
    }

    // the new prime field.
    inline void prime_changed() {
      par.precompute_tokens();
    }

    size_t n = 0;
  private:
    ShuntingYardParser par;
    //size_t integral_size, mi_size;
    //std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> amplitude_mapping;
  };

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
     *  @param vars_ A vector which specifies the variables of the functions. Each variable is allowed to be built of lower and upper case letters and numbers up to 16 characters.
     *  @param integral_families_ A vector which specifies the integral families of the amplitude.
     */
    AmplitudeParser(const std::vector<std::string>& vars_, const std::vector<std::string>& integral_families_);
    /**
     * Parses an expression from a file
     * @param amplitude_file the path to a file
     */
    void parse_amplitude_file(const std::string& amplitude_file);
    /**
     * Parses an expression from a string
     * @param amplitude_string an expression as a string
     */
    void parse_amplitude_string(std::string& amplitude_string);
    /**
     * Parses an expression from a file
     * @param amplitude_file the path to a file
     */
    void parse_file(const std::string& amplitude_file);
    /**
     * Parses an expression and splits it into factors and functions
     * @param amplitude an expression as a string
     * @return a vector of pairs of strings where the first entry is the function and the second the factor
     */
    std::vector<std::pair<std::string, std::string>> parse_string(std::string& amplitude, bool replacement = false);
    /**
     * Parses a replacement list from a file
     * @param ibp_table the path to a replacement list
     */
    void parse_ibp_table_file(const std::string& ibp_table);
    /**
     * Parses a replacement list from a string
     * @param ibp_table a replacement list as a string
     */
    void parse_ibp_table_string(std::string& ibp_table);
    /**
     * Checks for functions that do not have replacement rules
     * @return the number of functions without replacement rules
     */
    size_t check_for_unreplaced_masters();
    /**
     * Builds a black box for a given basis function
     * @param master the counter to indentify the basis function
     * @return an FFAmplitudeBlackBox object
     */
    FFAmplitudeBlackBox build_black_box(size_t master) const;
    /**
     * Returns an unsimplified coefficient of a basis function
     * @param master the counter to identify the basis function
     * @return the unsimplified coefficient of a basis function
     */
    std::string get_unsimplified_coef(size_t master) const;
    /**
     * Returns the ith basis function
     * @param i the counter to identify the basis function
     * @return the ith basis function
     */
    std::string get_master(size_t i);
  private:
    std::vector<std::string> vars; /**< Vector of variables*/
    std::vector<std::string> integral_families{}; /**< Vector of functions*/
    std::unordered_map<std::string, size_t> integrals{}; /**< Map of identified functions*/
    std::unordered_set<std::string> repl_integrals{}; /**< Map of functions that get replaced*/
    std::unordered_map<std::string, size_t> masters{}; /**< Map of basis functions*/
    std::unordered_map<size_t, std::string> masters_inv{}; /**< Inverse map of functions that get replaced*/
    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> amplitude_mapping{}; /**< Map that specifies how an expression is mapped via the found replacements*/
    std::vector<std::string> functions{}; /**< Vector of identified functions*/
    size_t distinct_integral_counter = 0; /**< A counter for distinct functions*/
    size_t distinct_master_counter = coef_type::COEF_TYPE_SIZE; /**< A counter for distinct basis functions*/
    size_t parser_counter = 0; /**< A counter to identify factors*/
    std::ofstream logger; /**< A logger file*/
  };
}
