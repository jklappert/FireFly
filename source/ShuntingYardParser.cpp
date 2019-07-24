#include "ShuntingYardParser.hpp"
#include "Logger.hpp"
#include <fstream>
#include <regex>
#include <chrono>

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
    auto time0 = std::chrono::high_resolution_clock::now();
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
    precompute_tokens();

    istream.close();
    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Parsed " + std::to_string(functions.size()) + " functions in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s.");
  }

  void ShuntingYardParser::parse(const std::string& fun_, bool use_regex) {
    std::string fun = fun_;

    if (use_regex) {
      for (const auto & el : vars_conv_map) {
        fun = std::regex_replace(fun, std::regex(el.second), std::string(1, el.first));
      }
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
      else if (is_variable(*l_ptr))
        tmp += *l_ptr;
      else if (is_operator(*l_ptr)) {
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

    parse(fun, false);

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
          // check then if number has more than 18 digits
          if (token.length() > 18)
            nums.push(FFInt(mpz_class(token)));
          else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);
              const char* var = tmp.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                nums.push(-values[vars_map.at(vars_conv_map.at(var[0]))]);
              else
                nums.push(-FFInt(std::stoull(tmp)));
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);
              const char* var = tmp.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                nums.push(values[vars_map.at(vars_conv_map.at(var[0]))]);
              else
                nums.push(FFInt(std::stoull(tmp)));
            } else {
              const char* var = token.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                nums.push(values[vars_map.at(vars_conv_map.at(var[0]))]);
              else
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

  std::vector<FFInt> ShuntingYardParser::evaluate_pre(const std::vector<FFInt>& values) const {
    std::vector<FFInt> res;
    res.reserve(functions.size());

    for (const auto & tokens : precomp_tokens) {
      std::stack<FFInt> nums;

      for (const auto & token : tokens) {
        switch (token.first) {
          case operands::OPERATOR : {
            // Pop two numbers
            FFInt a = nums.top();
            nums.pop();
            FFInt b = nums.top();
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
                nums.push(b.pow(a));
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
            nums.push(-values[token.second.n]);
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
        ERROR_MSG("Error in functional evaluation! Check your input.");
        std::exit(-1);
      }
    }

    return res;
  }

  int ShuntingYardParser::get_weight(const char c) const {
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

  bool ShuntingYardParser::is_operand(const char c) const {
    if (!is_operator(c) && c != '(' && c != ')' && vars_conv_map.find(c) == vars_conv_map.end() && c != ' ')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_operator(const char c) const {
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_variable(const char c) const {
    if (vars_conv_map.find(c) != vars_conv_map.end())
      return true;

    return false;
  }

  std::vector<std::vector<std::string>> ShuntingYardParser::get_rp_functions() const {
    return functions;
  }

  char ShuntingYardParser::get_var(int index) {
    if (int_var_map.empty())
      int_var_map = init_int_var_map();

    return int_var_map.at(index);
  }

  bool ShuntingYardParser::empty() const {
    return functions.empty();
  }

  void ShuntingYardParser::precompute_tokens() {
    precomp_tokens.clear();
    uint64_t size = functions.size();
    precomp_tokens = std::vector<std::vector<std::pair<uint8_t, FFInt>>> (size);

    for (uint64_t i = 0; i < size; ++i) {
      std::vector<std::string> tokens = functions[i];
      uint64_t t_size = tokens.size();
      precomp_tokens[i] = std::vector<std::pair<uint8_t, FFInt>> (t_size);

      for (uint64_t j = 0; j < t_size; ++j) {
        std::string token = tokens[j];

        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^") {
          switch (token[0]) {
            case '+': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::PLUS};
              break;
            }

            case '-': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::MINUS};
              break;
            }

            case '*': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::MULT};
              break;
            }

            case '/': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::DIV};
              break;
            }

            case '^': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::POW};
              break;
            }
          }
        } else {
          // check then if number has more than 18 digits
          if (token.length() > 18)
            precomp_tokens[i][j] = {operands::NUMBER, FFInt(mpz_class(token))};
          else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);
              const char* var = tmp.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                precomp_tokens[i][j] = {operands::NEG_VARIABLE, vars_map.at(vars_conv_map.at(var[0]))};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (-FFInt(std::stoull(tmp)))};
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);
              const char* var = tmp.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map.at(vars_conv_map.at(var[0]))};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(tmp)))};
            } else {
              const char* var = token.c_str();

              if (vars_conv_map.find(var[0]) != vars_conv_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map.at(vars_conv_map.at(var[0]))};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(token)))};
            }
          }
        }
      }
    }
  }
}
