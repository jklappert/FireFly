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

#include "firefly/config.hpp"
#include "firefly/FFInt.hpp"
#include "firefly/Logger.hpp"
#include "firefly/UintHasher.hpp"

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
     *  @param fun_ The rational function to be parsed as a string
     *  @param vars A vector which specifies the variables of the function
     *  @param validate_fun validates the function if set to true
     */
    void parse_function(std::string& fun_, const std::vector<std::string>& vars, bool validate_fun = false);
    /**
     *  TODO
     *  @param fun The rational function to be parsed as a string
     *  @param no_duplicates If true, do not add the function if it is already in the list, requires check_is_equal = true
     *  @return The position of the function in the list
     */
    size_t add_otf(const std::string & fun, const bool no_duplicates = false);
    /**
     *  Adds a function in RPN to the parser
     *  @param rpn_fun the function in RPN
     *  @return The position of the function in the list
     */
    size_t add_otf(const std::vector<std::string>& rpn_fun);
    /**
     *  Reserves memory
     *  @param number_of_functions how many entries shall be reserved
     */
    void reserve(size_t number_of_functions);
    /**
     *  Evaluates all functions for a given parameter point and returns their result.
     *  @param values A vector of FFIntTemp objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector.
     */
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate(const std::vector<FFIntTemp>& values) const;
    /**
     *  Evaluates all functions for a given parameter point and returs their result using precomputed values in bunches. This is in general much faster than ShuntingYardParser::evaluate.
     *  @param values A vector of vectors of FFIntTemp objects at which the parsed functions should be evaluated.
     *  @return The values of the parsed functions as a vector of vectors.
     */
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate_pre(const std::vector<FFIntTemp>& values) const noexcept;
    template<typename FFIntTemp>
    std::vector<FFIntTemp> evaluate_pre_2(const std::vector<FFIntTemp>& values) const noexcept;
    /**
     *  Returns the reverse polish notation of the parsed functions
     *  @return A vector of all parsed functions in reverse polish notation. Note the additional operators '~', '!', and ';'.
     */
    std::vector<std::vector<std::string>> get_rp_functions() const;
    /**
     *  Returns the reverse polish notation of the ith parsed function
     *  @param i the counter for the corresponding function
     *  @return A the ith function.
     */
    std::vector<std::string> get_rp_function(size_t i) const;
    /**
     *  Returns references the reverse polish notation of the parsed functions
     *  @return A vector of all parsed function references.
     */
    const std::vector<std::vector<std::string>> *get_rp_functions_ref() const;
    /**
     *  Checks if functions are stored in this class
     *  @return True if no functions are stored in this class
     */
    bool empty() const;
    /**
     *  Precomputes the tokes over the current prime field to be more efficient in evaluations
     *  @param froce enforces a precompute
     */
    void precompute_tokens(bool force = false);
    /**
     *  TODO
     *  @param elements_to_keep The elements which should be kept
     *  @return A map indicating on which positions the kept elements are know
     */
    std::unordered_map<size_t, size_t> trim(const std::unordered_set<size_t> & elements_to_keep);
  private:
    std::vector<std::vector<std::string>> functions {}; /**< This vector holds the input function in reverse polish notation where each operand/operator is separated */
    std::unordered_map<std::string, int> vars_map {}; /**< This map holds the conversion of the used variables to an integer value to map variables to parameter points */
    std::vector<std::vector<std::pair<uint8_t, FFInt>>> precomp_tokens {}; /**< This vector holds a collection of precomputed tokens where each operand/operator as a string is already converted to an FFInt */
    bool check_is_equal = false; /**< Indicates that functions are checked whether they are equal. Modifies the evaluation procedure */
    size_t prime_counter = 0; // TODO
    std::vector<FFInt> check_vars_1 {}; // TODO clear or keep?
    std::vector<FFInt> check_vars_2 {}; // TODO clear or keep?
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t, UintPairHasher> check_map {}; // TODO clear or keep?
    std::vector<uint64_t> evaluation_positions {}; /**< Stores the position of an evaluated function */
    std::uint64_t prime_internal = 0; /**< Stores the current prime and precomputes tokens only if this number has actually changed */
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
     *  @param fun The function which should be converted
     */
    std::vector<std::string> parse(std::string& fun);
    /**
     *  Checks expression for redundant parenthesis, removes them and throws an error if one encounters a mismatch of parenthesis
     *  @param r the expression as a string
     *  @param exp_n the expression number
     */
    void validate(std::string& r, uint64_t exp_n);
    /**
     *  Exits the program since a variable was not declared
     *  @param var the undeclared variable
     */
    void throw_not_declared_var_err(const std::string& var) const;
  };

  namespace tokens {
    enum tokens {
      PLUS,
      MINUS,
      MULT,
      DIV,
      POW,
      NEG_POW,
      POW_NEG,
      NEG_POW_NEG,
      OPERATOR,
      VARIABLE,
      NEG_VARIABLE,
      NUMBER
    }; /**< The allowed operands as an enum for fast identification. Note that ! means x ! 2 = -x^2 (! == POW_NEG) and ~ identifies a negative exponent (~ == NEG_POW), x ; 2 = -x^(-2) (; == NEG_POW_NEG) */
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

          // Evaluate and push the result back to the stack
          switch (token[0]) {
            case '+': {
              nums.top() += a;
              break;
            }

            case '-': {
              nums.top() -= a;
              break;
            }

            case '*': {
              nums.top() *= a;
              break;
            }

            case '/': {
              nums.top() /= a;
              break;
            }

            case '^': {
              nums.top() = std::move(nums.top().pow(a));
              break;
            }

            case '!': {
              nums.top() = std::move(-nums.top().pow(a));
              break;
            }

            case '~': {
              nums.top() = std::move(nums.top().pow(a.to_neg_int()));
              break;
            }

            case ';': {
              nums.top() = std::move(-nums.top().pow(a.to_neg_int()));
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
  std::vector<FFIntTemp> ShuntingYardParser::evaluate_pre(const std::vector<FFIntTemp>& values) const noexcept {
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
          case tokens::PLUS: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() += a;
            break;
          }

          case tokens::MINUS: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() -= a;
            break;
          }

          case tokens::MULT: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() *= a;
            break;
          }

          case tokens::DIV: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() /= a;
            break;
          }

          case tokens::POW: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() = std::move(pow(nums.top(), a));
            break;
          }

          case tokens::NEG_POW: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() = std::move(nums.top().pow(a.to_neg_int()));
            break;
          }

          case tokens::POW_NEG: {
            FFIntTemp a = nums.top();
            nums.pop();
            nums.top() = std::move(-pow(nums.top(), a));
            break;
          }

          case tokens::NEG_POW_NEG: {
            FFIntTemp a = nums.top();
            nums.pop();
	    nums.top() = std::move(-nums.top().pow(a.to_neg_int()));
            break;
          }

          case tokens::VARIABLE: {
            nums.push(values[token.second.n]);
            break;
          }

          case tokens::NEG_VARIABLE: {
            nums.push(neg_values[token.second.n]);
            break;
          }

          case tokens::NUMBER: {
            nums.push(token.second);
	    break;
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
  std::vector<FFIntTemp> ShuntingYardParser::evaluate_pre_2(const std::vector<FFIntTemp>& values) const noexcept {
    std::vector<FFIntTemp> res;
    res.reserve(functions.size());
    std::vector<FFIntTemp> neg_values;
    neg_values.reserve(values.size());

    for (const auto& el : values) {
      neg_values.emplace_back(-el);
    }

    for (const auto& tokens : precomp_tokens) {
      FFIntTemp* stack = new FFIntTemp[tokens.size()];//TODO evaluate maximum stack size. using the heap will make the evaluation slower!
      size_t stack_depth = 0;

      for (const auto& token : tokens) {
        switch (token.first) {
          case tokens::PLUS: {
	    stack[stack_depth - 1] += stack[stack_depth];
	    --stack_depth;
            break;
          }

          case tokens::MINUS: {
	    stack[stack_depth - 1] -= stack[stack_depth];
	    --stack_depth;
            break;
          }

          case tokens::MULT: {
	    stack[stack_depth - 1] *= stack[stack_depth];
	    --stack_depth;
            break;
          }

          case tokens::DIV: {
	    stack[stack_depth - 1] /= stack[stack_depth];
	    --stack_depth;
            break;
          }

          case tokens::POW: {
	    stack[stack_depth - 1] = std::move(pow(stack[stack_depth - 1], stack[stack_depth]));
	    --stack_depth;
            break;
          }

          case tokens::NEG_POW: {
	    stack[stack_depth - 1] = std::move(stack[stack_depth - 1].pow(stack[stack_depth].to_neg_int()));
	    --stack_depth;
            break;
          }

          case tokens::POW_NEG: {
	    stack[stack_depth - 1] = std::move(-pow(stack[stack_depth - 1], stack[stack_depth]));
	    --stack_depth;
            break;
          }

          case tokens::NEG_POW_NEG: {
	    stack[stack_depth - 1] = std::move(-stack[stack_depth - 1].pow(stack[stack_depth].to_neg_int()));
	    --stack_depth;
            break;
          }

          case tokens::VARIABLE: {
	    stack[++stack_depth] = values[token.second.n];
            break;
          }

          case tokens::NEG_VARIABLE: {
	    stack[++stack_depth] = neg_values[token.second.n];
            break;
          }

          case tokens::NUMBER: {
	    stack[++stack_depth] = token.second;
	    break;
          }
        }
      }

      if (stack_depth == 1)
        res.emplace_back(stack[1]);
      else {
	delete[] stack;
        ERROR_MSG("Error in functional evaluation! Please check your input.");
        std::exit(EXIT_FAILURE);
      }
      delete[] stack;
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
