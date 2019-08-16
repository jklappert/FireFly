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

#include <unordered_map>
#include <unordered_set>

namespace firefly {
  /**
   * @class ShuntingYardParser
   * @brief A parser for converting rational functions to reverse polish notation and evaluating the results
   */
  class ShuntingYardParser {
  public:
    /**
     *  Constructor which parses a list of rational functions and prepares them for evaluation.
     *  @param file The path to the file in which the rational functions are stored
     *  @param vars A vector which specifies the variables of the functions. Each variable is allowed to be built of lower and upper case letters and numbers up to 16 characters.
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
     *  Parses a rational function given as a string for internal functions only
     *  @param fun The rational function to be parsed as a string
     *  @param vars A vector which specifies the variables of the function.
     */
    void parse_function_internal(const std::string& fun, const std::vector<std::string>& vars);
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
     *  Checks if functions are stored in this class
     *  @return True if no functions are stored in this class
     */
    bool empty() const;
    /**
     *  Precomputes the tokes over the current prime field to be more efficient in evaluations
     */
    void precompute_tokens();
  private:
    std::vector<std::vector<std::string>> functions {}; /**< This vector holds the input function in reverse polish notation where each operand/operator is separated */
    std::unordered_map<std::string, int> vars_map {}; /**< This map holds the conversion of the used variables to an integer value to map variables to parameter points */
    std::vector<std::vector<std::pair<uint8_t, FFInt>>> precomp_tokens {}; /**< This vector holds a collection of precomputed tokens where each operand/operator as a string is already converted to an FFInt */
    /**
     *  Initializes the conversion map
     *  @return The conversion map
     */
    static const std::unordered_set<char> chars;/**< The collection of supported characters for variable definition */
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
     *  Converts a function in reverse polish notation for internal functions
     *  @param fun_ The function which should be converted
     */
    void parse_internal(const std::string& fun_);
    /**
     *  Checks expression for redundant parenthesis, removes them and throws an error if one encounters a mismatch of parenthesis
     *  @param line the expression as a string
     *  @param exp_n the expression number
     */
    std::string validate(const std::string& line, uint32_t exp_n);
    /**
     *  Exits the program since a variable was not declared
     *  @param var the undeclared variable
     */
    void throw_not_declared_var_err(const std::string& var) const;
  };

  namespace operators {
    enum operators {
      PLUS,
      MINUS,
      MULT,
      DIV,
      POW
    }; /**< The allowed operators as an enum for fast identification */
  }

  namespace operands {
    enum operands {
      OPERATOR,
      VARIABLE,
      NEG_VARIABLE,
      NUMBER
    }; /**< The allowed operands as an enum for fast identification */
  }
}
