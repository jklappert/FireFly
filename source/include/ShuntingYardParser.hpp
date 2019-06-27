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

namespace firefly{

  class ShuntingYardParser {
  public:
    /**
     *  Constructor which parses a list of rational functions and prepares them for evaluation.
     *  @param file The path to the file in which the rational functions are stored
     *  @param vars A vector which specifies the variables of the functions.
     */
    ShuntingYardParser(std::string file, std::vector<std::string> vars);
    /**
     *  Evaluates all functions for a given parameter point and returs their result.
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    std::vector<FFInt> evaluate(const std::vector<FFInt>& values);
    /**
     *  Returns the reverse polish notation of the parsed functions
     *  @return A vector of all parsed functions in reverse polish notation
     */
    std::vector<std::vector<std::string>> get_rp_functions();
  private:
    std::vector<std::vector<std::string>> functions {};
    std::unordered_map<std::string, int> vars_map {};
    std::unordered_map<char, std::string> vars_conv_map {};
    std::unordered_map<int, char> int_var_map  = {{0,'a'},{1,'b'},{2,'c'},{3,'d'},
    {4,'e'},{5,'f'},{6,'g'},{7,'h'},{8,'i'},{9,'j'},{10,'k'},{11,'l'},{12,'m'},
    {13,'n'},{14,'o'},{15,'p'},{16,'q'},{17,'r'},{18,'s'},{19,'t'},{20,'u'},{21,'v'},
    {22,'w'},{23,'x'},{24,'y'},{25,'z'}};
    int get_weight(const char c);
    bool is_operator(char c);
    bool is_operand(char c);
    bool is_variable(char c);
  };
}