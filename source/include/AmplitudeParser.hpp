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

#include "BlackBoxBase.hpp"
#include "FFInt.hpp"
#include "Logger.hpp"
#include "ShuntingYardParser.hpp"

#include <unordered_map>
#include <unordered_set>
#include <list>

namespace firefly {
  namespace coef_type {
    enum coef_type {
      PREFACTOR = 0
    }; /**< The key to identify a prefactor from a MI coefficient */
  }

  /**
   * @class FFAmplitudeBlackBox
   * @brief The black box of the amplitude parser which inserts IBP relations into an amplitude and returns the reduced amplitude in terms of master integrals
   */
  class FFAmplitudeBlackBox : public BlackBoxBase<FFAmplitudeBlackBox> {
  public:
    FFAmplitudeBlackBox(const size_t integral_size_,
                        const size_t mi_size_,
                        const std::vector<std::string>& vars,
                        const std::vector<std::string>& functions,
                        const std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>& amplitude_mapping_) {
      integral_size = integral_size_;
      mi_size = mi_size_ - 1;
      par = ShuntingYardParser(functions, vars, true);
      n = vars.size();
      amplitude_mapping = amplitude_mapping_;
    };

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      auto res = par.evaluate_pre(values);

      std::vector<FFIntTemp> result (mi_size, FFIntTemp(0));

      for (size_t i = 0; i != integral_size; ++i) {
        FFIntTemp prefactor(0);
        for (const auto& mi_map : amplitude_mapping[i]) {
          bool pref_only = true;
          if (mi_map.second == coef_type::PREFACTOR)
            prefactor += res[mi_map.first];
          else {
            pref_only = false;
            result[mi_map.second - 1] += prefactor*res[mi_map.first];
          }
        }
      }

      return result;
    }

    // the new prime field.
    inline void prime_changed() {
      par.precompute_tokens();
    }

    size_t n = 0;
  private:
    ShuntingYardParser par;
    size_t integral_size, mi_size;
    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> amplitude_mapping;
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
     * 
     */
    void parse_amplitude_file(const std::string& amplitude_file);
    /**
     * 
     */
    void parse_amplitude_string(const std::string& amplitude_string);
    /**
     * 
     */
    void parse_file(const std::string& amplitude_file);
    /**
     * 
     */
    std::vector<std::pair<std::string, std::string>> parse_string(const std::string& amplitude);
    /**
     * 
     */
    void parse_ibp_table_file(const std::string& ibp_table);
    /**
     * 
     */
    void parse_ibp_table_string(const std::string& ibp_table);
    /**
     * 
     */
    FFAmplitudeBlackBox build_black_box();
    /**
     * 
     */
    std::string get_master(size_t i);
  private:
    std::vector<std::string> vars;
    std::vector<std::string> integral_families{};
    std::unordered_map<std::string, size_t> integrals{};
    std::unordered_set<std::string> repl_integrals{};
    std::unordered_map<std::string, size_t> masters{};
    std::unordered_map<size_t, std::string> masters_inv{};
    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> amplitude_mapping{};
    std::vector<std::string> functions{};
    ShuntingYardParser sy_parser;
    size_t distinct_integral_counter = 0;
    size_t distinct_master_counter = 1;
    size_t parser_counter = 0;
  };
}
