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
#include <unordered_set>

namespace firefly {
  /**
   * @class ShuntingYardParser
   * @brief A functional parser using reverse polish notation
   */
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
     *  Evaluates all functions for a given parameter point and returs their result using precomputed values. This is in general much faster than ShuntingYardParser::evaluate.
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    std::vector<FFInt> evaluate_pre(const std::vector<FFInt>& values) const;
    /**
     *  Evaluates all functions for a given parameter point and returs their result using precomputed values in bunches. This is in general much faster than ShuntingYardParser::evaluate.
     *  @param values A vector of vectors of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector of vectors.
     */
    std::vector<std::vector<FFInt>> evaluate_pre(const std::vector<std::vector<FFInt>>& values) const;
    /**
     *  Returns the reverse polish notation of the parsed functions
     *  @return A vector of all parsed functions in reverse polish notation with changed variable names according to int_var_map
     */
    std::vector<std::vector<std::string>> get_rp_functions() const;
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
    bool empty() const;
    /**
     *  Precomputes the tokes over the current prime field to be more efficient in evaluations
     */
    void precompute_tokens();
  private:
    std::vector<std::vector<std::string>> functions {};
    std::unordered_map<std::string, int> vars_map {};
    std::vector<std::vector<std::pair<uint8_t, FFInt>>> precomp_tokens {};
    /**
     *  Initializes the conversion map
     *  @return The conversion map
     */
    std::unordered_set<char> chars {{'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'}};
    /**
     *  Gets the weight of an operator
     *  @param c The operator as a character
     *  @return An integer which marks the weight
     */
    int get_weight(const char c) const;
    /**
     *  Checks if a character is an operator
     *  @param c The character which should be checked
     *  @return True if the charater is an operator
     */
    bool is_operator(const char c) const;
    /**
     *  Checks if a character is an operand
     *  @param c The character which should be checked
     *  @return True if the charater is an operand
     */
    bool is_operand(const char c) const;
    /**
     *  Checks if a character is a variable
     *  @param c The character which should be checked
     *  @return True if the charater is a variable
     */
    bool is_variable(const char c) const;
    /**
     *  Converts a function in reverse polish notation
     *  @param fun_ The function which should be converted
     */
    void parse(const std::string& fun_);
    /**
     *  Checks expression for redundant parenthesis, removes them and throws an error if one encounters a mismatch of parenthesis
     *  @param line the expression as a string
     *  @param exp_n the expression number
     */
    std::string validate(const std::string& line, uint32_t exp_n);
  };

  namespace operators {
    enum operators {
      PLUS,
      MINUS,
      MULT,
      DIV,
      POW
    };
  }

  namespace operands {
    enum operands {
      OPERATOR,
      VARIABLE,
      NEG_VARIABLE,
      NUMBER
    };
  }
}