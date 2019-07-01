#include "ShuntingYardParser.hpp"
#include "Logger.hpp"
#include <fstream>
#include <regex>

namespace firefly {

  std::unordered_map<int, char> ShuntingYardParser::int_var_map;

  ShuntingYardParser::ShuntingYardParser() {
    if (int_var_map.empty())
      int_var_map = init_int_var_map();
  }

  ShuntingYardParser::ShuntingYardParser(const std::string& file, const std::vector<std::string>& vars) {

    if (int_var_map.empty())
      int_var_map = init_int_var_map();

    INFO_MSG("Parsing functions in '" + file + "'.");
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
      vars_conv_map.emplace(std::make_pair(int_var_map.at(i), vars[i]));
    }

    while (std::getline(istream, line)) {
      parse(line);
    }

    functions.shrink_to_fit();

    istream.close();
    INFO_MSG("Parsed " + std::to_string(functions.size()) + " functions.");
  }

  void ShuntingYardParser::parse(const std::string& fun_) {
    std::string fun = fun_;

    for (const auto & el : vars_conv_map) {
      fun = std::regex_replace(fun, std::regex(el.second), std::string(1, el.first));
    }

    char const* l_ptr = fun.c_str();
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

        if (!op_stack.empty() && *(l_ptr - 1) == '(') {
          tmp.insert(tmp.begin(), *l_ptr);
        } else if (op_stack.empty() && pf.empty()) {
          tmp.insert(tmp.begin(), *l_ptr);
        } else {

          while (!op_stack.empty() && op_stack.top() != '(' && get_weight(op_stack.top()) >= get_weight(*l_ptr)) {
            pf.emplace_back(std::string(1, op_stack.top()));
            op_stack.pop();
          }

          op_stack.push(*l_ptr);
        }
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

  void ShuntingYardParser::parse_function(const std::string& fun, const std::vector<std::string>& vars) {
    if (int_var_map.empty())
      int_var_map = init_int_var_map();

    for (uint32_t i = 0; i < vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
      vars_conv_map.emplace(std::make_pair(int_var_map.at(i), vars[i]));
    }

    parse(fun);

    functions.shrink_to_fit();
  }

  std::vector<FFInt> ShuntingYardParser::evaluate(const std::vector<FFInt>& values) const {
    std::vector<FFInt> res;
    res.reserve(functions.size());

    for (const auto & tokens : functions) {
      std::stack<FFInt> nums;

      for (const auto & token : tokens) {
        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^") {
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
        } else {
          const char* var = token.c_str();

          // check first if the token is a is_variable
          if (vars_conv_map.find(var[0]) != vars_conv_map.end())
            nums.push(values[vars_map.at(vars_conv_map.at(var[0]))]);
          // check then if number has more than 18 digits
          else if (token.length() > 18)
            nums.push(FFInt(mpz_class(token)));
          else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);
              nums.push(-FFInt(std::stoull(tmp)));
            } else {
              nums.push(FFInt(std::stoull(token)));
            }
          }
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

  std::vector<std::vector<std::string>> ShuntingYardParser::get_rp_functions() {
    return functions;
  }

  char ShuntingYardParser::get_var(int index) {
    if (int_var_map.empty())
      int_var_map = init_int_var_map();

    return int_var_map.at(index);
  }

  bool ShuntingYardParser::empty() {
    return functions.empty();
  }

  void precomp_token(const std::string& token) {
    /*if (token.length() > 18)
      precomp.emplace(FFInt(mpz_class(token)))
    else {
      if (token[0] == '-') {
        std::string tmp = token;
        tmp.erase(0, 1);
        nums.push(-FFInt(std::stoull(tmp)));
      } else {
        nums.push(FFInt(std::stoull(token)));
      }*/
    }
  }
