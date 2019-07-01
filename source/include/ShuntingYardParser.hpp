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
#include <stack>
#include <unordered_map>

namespace firefly {

  class ShuntingYardParser {
  public:
    /**
     *  Constructor which parses a list of rational functions and prepares them for evaluation.
     *  @param file The path to the file in which the rational functions are stored
     *  @param vars A vector which specifies the variables of the functions.
     */
    ShuntingYardParser(const std::string& file, const std::vector<std::string>& vars);
    /**
     *  Default constructor
     */
    ShuntingYardParser();
    /**
     *  Parses a rational function given as a string
     *  @param fun The rational function to be parsed as a string
     *  @param vars A vector which specifies the variables of the function.
     */
    void parse_function(const std::string& fun, const std::vector<std::string>& vars);
    /**
     *  Evaluates all functions for a given parameter point and returs their result.
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    std::vector<FFInt> evaluate(const std::vector<FFInt>& values) const;
    /**
     *  Returns the reverse polish notation of the parsed functions
     *  @return A vector of all parsed functions in reverse polish notation with changed variable names according to int_var_map
     */
    std::vector<std::vector<std::string>> get_rp_functions();
    /**
     *  Returns the mapped variable to a given index
     *  @param index The index of the variable
     *  @return The mapped variable
     */
    static char get_var(int index);
    /**
     *  Checks if functions are stored in this class
     *  @return True if no functions are stored in this class
     */
    bool empty();
  private:
    std::unordered_map<std::string, FFInt> precomp {};
    std::vector<std::vector<std::string>> functions {};
    std::unordered_map<std::string, int> vars_map {};
    std::unordered_map<char, std::string> vars_conv_map {};
    void precomp_token(const std::string& token);
    static std::unordered_map<int, char> int_var_map;
    /**
     *  Initializes the conversion map
     *  @return The conversion map
     */
    static std::unordered_map<int, char> init_int_var_map() {
      std::unordered_map<int, char> m = {{0, 'a'}, {1, 'b'}, {2, 'c'}, {3, 'd'},
        {4, 'e'}, {5, 'f'}, {6, 'g'}, {7, 'h'}, {8, 'i'}, {9, 'j'}, {10, 'k'}, {11, 'l'}, {12, 'm'},
        {13, 'n'}, {14, 'o'}, {15, 'p'}, {16, 'q'}, {17, 'r'}, {18, 's'}, {19, 't'}, {20, 'u'}, {21, 'v'},
        {22, 'w'}, {23, 'x'}, {24, 'y'}, {25, 'z'}
      };
      return m;
    }
    /**
     *  Gets the weight of an operator
     *  @param c The operator as a character
     *  @return An integer which marks the weight
     */
    int get_weight(const char c);
    /**
     *  Checks if a character is an operator
     *  @param c The character which should be checked
     *  @return True if the charater is an operator
     */
    bool is_operator(char c);
    /**
     *  Checks if a character is an operand
     *  @param c The character which should be checked
     *  @return True if the charater is an operand
     */
    bool is_operand(char c);
    /**
     *  Checks if a character is a variable
     *  @param c The character which should be checked
     *  @return True if the charater is a variable
     */
    bool is_variable(char c);
    /**
     *  Converts a function in reverse polish notation
     *  @param fun The function which should be converted
     */
    void parse(const std::string&);
  };
}
