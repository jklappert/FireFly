#include "ShuntingYardParser.hpp"
#include "Logger.hpp"
#include <fstream>
#include <regex>

namespace firefly {

  ShuntingYardParser::ShuntingYardParser(std::string file, std::vector<std::string> vars) {

    // Check if file exists
    std::ifstream infile(file);

    if (!infile.good()) {
      ERROR_MSG("File '" + file + "' does not exist!");
      std::exit(-1);
    }

    std::ifstream istream;
    istream.open(file);

    std::string line;

    for (uint32_t i = 0; i < vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
      vars_conv_map.emplace(std::make_pair(int_var_map[i], vars[i]));
    }

    while (std::getline(istream, line)) {
      for (const auto & el : vars_conv_map) {
        line = std::regex_replace(line, std::regex(el.second), std::string(1, el.first));
      }

      char const* l_ptr = line.c_str();
      std::string tmp = ""; // Used for numbers
      std::vector<std::string> pf = {};
      std::stack<char> op_stack;

      // Pick one character at a time until we reach the end of the line
      while (*l_ptr != '\0') {
        // If operand, add it to postfix string
        // If operator pop operators off the stack until it is empty
        if (is_operand(*l_ptr))
          tmp += *l_ptr;
        else if (is_variable(*l_ptr)) {
          if (tmp.length() > 0) {
            pf.emplace_back(tmp);
            tmp = "";
          }

          pf.emplace_back(std::string(1, *l_ptr));
        } else if (is_operator(*l_ptr)) {
          if (tmp.length() > 0) {
            pf.emplace_back(tmp);
            tmp = "";
          }

          while (!op_stack.empty() && op_stack.top() != '(' && get_weight(op_stack.top()) >= get_weight(*l_ptr)) {
            pf.emplace_back(std::string(1, op_stack.top()));
            op_stack.pop();
          }

          op_stack.push(*l_ptr);
        }
        // Push all open parenthesis to the stack
        else if (*l_ptr == '(')
          op_stack.push(*l_ptr);
        // When reaching a closing one, pop off operators from the stack until an opening one is found
        else if (*l_ptr == ')') {
          if (tmp.length() > 0) {
            pf.emplace_back(tmp);
            tmp = "";
          }

          while (!op_stack.empty()) {
            if (op_stack.top() == '(') {
              op_stack.pop();
              break;
            }

            pf.emplace_back(std::string(1, op_stack.top()));
            op_stack.pop();
          }
        }

        // Proceed to the next character
        l_ptr ++;
      }

      if (tmp.length() > 0)
        pf.emplace_back(tmp);

      while (!op_stack.empty()) {
        pf.emplace_back(std::string(1, op_stack.top()));
        op_stack.pop();
      }

      functions.emplace_back(pf);
    }

    functions.shrink_to_fit();

    istream.close();
  }

  std::vector<FFInt> ShuntingYardParser::evaluate(const std::vector<FFInt>& values) {
    std::vector<FFInt> res;
    res.reserve(functions.size());

    for (const auto & tokens : functions) {
      std::stack<FFInt> nums;

      for (const auto & token : tokens) {
        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^") {
          // Usual case where two numbers are connected by an operator
          if (nums.size() > 1) {
            // Pop two numbers
            FFInt a = nums.top();
            nums.pop();
            FFInt b = nums.top();
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
            }
          }
          // Needed for negative numbers where one number is connected with an operator
          else {
            // Pop a number
            FFInt a = nums.top();
            nums.pop();

            // Evaluate and push the result back to the stack
            switch (token[0]) {
              case '+': {
                nums.push(a);
                break;
              }

              case '-': {
                nums.push(-a);
                break;
              }
            }
          }
        } else {
          const char* var = token.c_str();

          // check first if the token is a is_variable
          if (vars_conv_map.find(var[0]) != vars_conv_map.end())
            nums.push(values[vars_map[vars_conv_map[var[0]]]]);
          // check then if number has more than 18 digits
          else if (token.length() > 18)
            nums.push(FFInt(mpz_class(token)));
          else
            nums.push(FFInt(std::stoull(token)));
        }
      }

      if (nums.size())
        res.emplace_back(nums.top());
      else {
        ERROR_MSG("Error in functional evaluation! Check your input.");
        std::exit(-1);
      }
    }

    return res;
  }

  int ShuntingYardParser::get_weight(const char c) {
    switch (c) {
      case '^':
        return 3;

      case '/':
      case '*':
        return 2;

      case '+':
      case '-':
        return 1;

      default :
        return 0;
    }
  }

  bool ShuntingYardParser::is_operand(char c) {
    if (!is_operator(c) && c != '(' && c != ')' && vars_conv_map.find(c) == vars_conv_map.end() && c != ' ')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_operator(char c) {
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_variable(char c) {
    if (vars_conv_map.find(c) != vars_conv_map.end())
      return true;

    return false;
  }
}
