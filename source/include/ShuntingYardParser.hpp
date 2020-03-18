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
#include <stack>

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
     *  @param check_is_equal Sets the option to check whether parsed functions are equal to only evaluate them once
     */
    ShuntingYardParser(const std::string& file, const std::vector<std::string>& vars, bool check_is_equal = false);
    /**
     *  Constructor which parses a list of rational functions and prepares them for evaluation.
     *  @param funs A collection of functions
     *  @param vars A vector which specifies the variables of the functions. Each variable is allowed to be built of lower and upper case letters and numbers up to 16 characters.
     *  @param check_is_equal Sets the option to check whether parsed functions are equal to only evaluate them once
     */
    ShuntingYardParser(const std::vector<std::string>& funs, const std::vector<std::string>& vars, bool check_is_equal = false);
    /**
     *  Default constructor
     */
    ShuntingYardParser();
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
    // TODO
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate(const std::vector<FFIntTemp>& values) const;
    /**
     *  Evaluates all functions for a given parameter point and returns their result using precomputed values. This is in general much faster than ShuntingYardParser::evaluate.
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    /**
     *  Evaluates all functions for a given parameter point and returs their result using precomputed values in bunches. This is in general much faster than ShuntingYardParser::evaluate.
     *  @param values A vector of vectors of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector of vectors.
     */
    // TODO
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate_pre(const std::vector<FFIntTemp>& values) const;
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
    bool check_is_equal = false; /**< Indicates that functions are checked whether they are equal. Modifies the evaluation procedure */
    std::vector<uint64_t> evaluation_positions; /**< Stores the position of an evaluated function */
    /**
     *  Evaluates a single function for a given parameter point and returns their result.
     *  @param fun A function in reverse polish notation
     *  @param values A vector of FFInt objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    FFInt evaluate(const std::vector<std::string>& fun, const std::vector<FFInt>& values) const;
    /**
     *  Parses either a collection of strings or a file
     *  @param funs the collection of functions
     *  @param bool is_file indicates wether the collection is stored in a file
     *
     */
    void parse_collection(const std::vector<std::string>& funs, bool is_file);
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
     *  Checks expression for redundant parenthesis, removes them and throws an error if one encounters a mismatch of parenthesis
     *  @param line the expression as a string
     *  @param exp_n the expression number
     */
    std::string validate(const std::string& line, uint64_t exp_n);
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
      POW,
      POW_NEG
    }; /**< The allowed operators as an enum for fast identification. Note that ! means x ! 2 = -x^2 (! == POW_NEG) */
  }

  namespace operands {
    enum operands {
      OPERATOR,
      VARIABLE,
      NEG_VARIABLE,
      NUMBER
    }; /**< The allowed operands as an enum for fast identification */
  }

  template<typename FFIntTemp>
  std::vector<FFIntTemp> ShuntingYardParser::evaluate(const std::vector<FFIntTemp>& values) const {
    std::vector<FFIntTemp> res;
    res.reserve(functions.size());

    for (const auto& tokens : functions) {
      std::stack<FFIntTemp> nums;

      for (const auto& token : tokens) {
        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^" || token == "!") {
          // Pop two numbers
          FFIntTemp a = nums.top();
          nums.pop();
          FFIntTemp b = nums.top();
          nums.pop();

          // Evaluate and push the result back to the stack
          switch (token[0]) {
            case '+': {
              nums.push(a + b);
              break;
            }

            case '-': {
              nums.push(b - a);
              break;
            }

            case '*': {
              nums.push(b * a);
              break;
            }

            case '/': {
              nums.push(b / a);
              break;
            }

            case '^': {
              nums.push(b.pow(a));
              break;
            }

            case '!': {
              nums.push(-b.pow(a));
              break;
            }
          }
        } else {
          // check then if number has more than 18 digits
          if (token.length() > 18) {
            std::string tmp = token;

            if (token[0] == '+')
              tmp.erase(0, 1);

            nums.push(FFIntTemp(mpz_class(tmp)));
          } else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                nums.push(-values[vars_map.at(tmp)]);
              else if (isdigit(tmp[0]))
                nums.push(-FFIntTemp(std::stoull(tmp)));
              else
                throw_not_declared_var_err(tmp);
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                nums.push(values[vars_map.at(tmp)]);
              else  if (isdigit(tmp[0]))
                nums.push(FFIntTemp(std::stoull(tmp)));
              else
                throw_not_declared_var_err(tmp);
            } else {
              if (vars_map.find(token) != vars_map.end())
                nums.push(values[vars_map.at(token)]);
              else if (isdigit(token[0]))
                nums.push(FFIntTemp(std::stoull(token)));
              else
                throw_not_declared_var_err(token);
            }
          }
        }
      }

      if (nums.size())
        res.emplace_back(nums.top());
      else {
        ERROR_MSG("Error in functional evaluation! Please check your input.");
        std::exit(EXIT_FAILURE);
      }
    }

    if (!check_is_equal)
      return res;
    else {
      uint64_t s = evaluation_positions.size();

      if (functions.size() != s) {
        std::vector<FFIntTemp> full_res;
        full_res.reserve(s);

        for (const auto& el : evaluation_positions) {
          full_res.emplace_back(res[el]);
        }

        return full_res;
      } else
        return res;
    }
  }

  template<typename FFIntTemp>
  std::vector<FFIntTemp> ShuntingYardParser::evaluate_pre(const std::vector<FFIntTemp>& values) const {
    std::vector<FFIntTemp> res;
    res.reserve(functions.size());
    std::vector<FFIntTemp> neg_values;
    neg_values.reserve(values.size());

    for (const auto& el : values) {
      neg_values.emplace_back(-el);
    }

    for (const auto& tokens : precomp_tokens) {
      std::stack<FFIntTemp> nums;

      for (const auto& token : tokens) {
        switch (token.first) {
          case operands::OPERATOR : {
            // Pop two numbers
            FFIntTemp a = nums.top();
            nums.pop();
            FFIntTemp b = nums.top();
            nums.pop();

            switch (token.second.n) {
              case operators::PLUS: {
                nums.push(a + b);
                break;
              }

              case operators::MINUS: {
                nums.push(b - a);
                break;
              }

              case operators::MULT: {
                nums.push(b * a);
                break;
              }

              case operators::DIV: {
                nums.push(b / a);
                break;
              }

              case operators::POW: {
                nums.push(pow(b, a));
                break;
              }

              case operators::POW_NEG: {
                nums.push(-pow(b, a));
                break;
              }
            }

            break;
          }

          case operands::VARIABLE : {
            nums.push(values[token.second.n]);
            break;
          }

          case operands::NEG_VARIABLE : {
            //nums.push(-values[token.second.n]);
            nums.push(neg_values[token.second.n]);
            break;
          }

          case operands::NUMBER: {
            nums.push(token.second);
          }
        }
      }

      if (nums.size())
        res.emplace_back(nums.top());
      else {
        ERROR_MSG("Error in functional evaluation! Please check your input.");
        std::exit(EXIT_FAILURE);
      }
    }

    if (!check_is_equal)
      return res;
    else {
      uint64_t s = evaluation_positions.size();

      if (functions.size() != s) {
        std::vector<FFIntTemp> full_res;
        full_res.reserve(s);

        for (const auto& el : evaluation_positions) {
          full_res.emplace_back(res[el]);
        }

        return full_res;
      } else
        return res;
    }
  }
}
