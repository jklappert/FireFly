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

#include "ShuntingYardParser.hpp"
#include "BaseReconst.hpp"
#include "ReconstHelper.hpp"

#include <chrono>
#include <fstream>

namespace firefly {

  const std::unordered_set<char> ShuntingYardParser::chars = {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
    }
  };

  ShuntingYardParser::ShuntingYardParser() {}

  ShuntingYardParser::ShuntingYardParser(const std::string& file, const std::vector<std::string>& vars, bool check_is_equal_) {
    INFO_MSG("Parsing function(s) in '" + file + "'");
    check_is_equal = check_is_equal_;

    for (uint32_t i = 0; i != vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
    }

    std::vector<std::string> tmp = {file};
    parse_collection(tmp, true);
  }

  ShuntingYardParser::ShuntingYardParser(const std::vector<std::string>& funs, const std::vector<std::string>& vars, bool check_is_equal_) {
    INFO_MSG("Parsing collection of " + std::to_string(funs.size()) + " function(s)");
    check_is_equal = check_is_equal_;

    for (uint32_t i = 0; i != vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
    }

    parse_collection(funs, false);
  }

  void ShuntingYardParser::parse_collection(const std::vector<std::string>& funs, bool is_file) {
    size_t prime_counter = 0;
    std::vector<FFInt> check_vars_1;
    std::vector<FFInt> check_vars_2;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t, UintPairHasher> check_map;

    for (const auto& p : primes()) {
      if (FFInt::p == p)
        break;

      ++prime_counter;
    }

    if (check_is_equal) {
      size_t s = vars_map.size();
      check_vars_1.reserve(s);
      check_vars_2.reserve(s);
      BaseReconst base;
      uint64_t seed = static_cast<uint64_t>(std::time(0));
      base.set_seed(seed);

      FFInt::set_new_prime(primes()[prime_counter != 299 ? prime_counter + 1 : prime_counter - 1]);

      for (size_t i = 0; i != s; ++i) {
        check_vars_1.emplace_back(base.get_rand_64());
      }

      FFInt::set_new_prime(primes()[prime_counter]);

      for (size_t i = 0; i != s; ++i) {
        check_vars_2.emplace_back(base.get_rand_64());
      }
    }

    auto time0 = std::chrono::high_resolution_clock::now();
    uint64_t equal_fun_c = 0;
    uint64_t parsed_fun_c = 0;

    if (is_file) {
      std::string file = funs[0];
      // Check if file exists
      std::ifstream infile(file);

      if (!infile.good()) {
        ERROR_MSG("File '" + file + "' does not exist!");
        std::exit(EXIT_FAILURE);
      }

      std::ifstream istream;
      istream.open(file);

      std::string line;

      uint64_t line_c = 1;

      while (std::getline(istream, line, ';')) {
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());

        if (line.length() > 0) {
          line = validate(line, line_c);
          parse(line);
          ++parsed_fun_c;

          if (check_is_equal) {
            FFInt::set_new_prime(primes()[prime_counter != 299 ? prime_counter + 1 : prime_counter - 1]);
            FFInt v1 = evaluate(functions.back(), check_vars_1);
            FFInt::set_new_prime(primes()[prime_counter]);
            FFInt v2 = evaluate(functions.back(), check_vars_2);

            if (check_map.find(std::make_pair(v1.n, v2.n)) != check_map.end()) {
              functions.pop_back();
              ++equal_fun_c;
              evaluation_positions.emplace_back(check_map[std::make_pair(v1.n, v2.n)]);
            } else {
              uint64_t s = functions.size() - 1;
              check_map.emplace(std::make_pair(std::make_pair(v1.n, v2.n), s));
              evaluation_positions.emplace_back(s);
            }
          }
        }

	if (line_c + 1 < funs.size() + 1)
	  std::cerr << "\033[1;34mFireFly info:\033[0m " << line_c + 1 << " / " << funs.size() << "\r";
        line_c++;
      }

      istream.close();
    } else {
      for (size_t i = 0; i != funs.size(); ++i) {
        std::string fun = funs[i];
        fun = validate(fun, i);
        parse(fun);
        ++parsed_fun_c;

        if (check_is_equal) {
          FFInt::set_new_prime(primes()[prime_counter != 299 ? prime_counter + 1 : prime_counter - 1]);
          FFInt v1 = evaluate(functions.back(), check_vars_1);
          FFInt::set_new_prime(primes()[prime_counter]);
          FFInt v2 = evaluate(functions.back(), check_vars_2);

          if (check_map.find(std::make_pair(v1.n, v2.n)) != check_map.end()) {
            functions.pop_back();
            ++equal_fun_c;
            evaluation_positions.emplace_back(check_map[std::make_pair(v1.n, v2.n)]);
          } else {
            uint64_t s = functions.size() - 1;
            check_map.emplace(std::make_pair(std::make_pair(v1.n, v2.n), s));
            evaluation_positions.emplace_back(s);
          }
        }

	if (i + 1 < funs.size() + 1)
	  std::cerr << "\033[1;34mFireFly info:\033[0m " << i + 1 << " / " << funs.size() << "\r";
      }
    }

    functions.shrink_to_fit();
    precompute_tokens();

    auto time1 = std::chrono::high_resolution_clock::now();

    if (!check_is_equal) {
      INFO_MSG("Parsed " + std::to_string(parsed_fun_c) + " function(s) in "
               + std::to_string(std::chrono::duration<double>(time1 - time0).count())
               + " s");
    } else {
      evaluation_positions.shrink_to_fit();
      BaseReconst::reset();
      INFO_MSG("Parsed " + std::to_string(parsed_fun_c) + " function(s) in "
               + std::to_string(std::chrono::duration<double>(time1 - time0).count())
               + " s");
      INFO_MSG("Found " + std::to_string(parsed_fun_c - equal_fun_c) + " different function(s)");
    }
  }

  void ShuntingYardParser::parse(const std::string& fun_) {
    std::string fun = fun_;

    // Check for global signs
    if (fun.size() > 2 && (fun[0] == '+' || fun[0] == '-') && fun[1] == '(') {
      fun.insert(fun.begin(), '0');
    }

    char const* l_ptr = fun.c_str();
    std::string tmp = ""; // Used for numbers
    std::vector<std::string> pf = {};
    std::stack<char> op_stack;
    bool neg_pow = false;
    uint32_t counter = 0;

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

        // Check for cases like +(-(x+...))
        if (!op_stack.empty() && *(l_ptr - 1) == '(') {
          if (*(l_ptr + 1) == '(') {
            pf.emplace_back("0");
            op_stack.push(*l_ptr);
          } else
            tmp.insert(tmp.begin(), *l_ptr);
        } else if (op_stack.empty() && pf.empty())
          tmp.insert(tmp.begin(), *l_ptr);
        else {

          while (!op_stack.empty() && op_stack.top() != '(' && get_weight(op_stack.top()) >= get_weight(*l_ptr)) {
            pf.emplace_back(std::string(1, op_stack.top()));
            op_stack.pop();
          }

          // Check for negative exponents
          if (*l_ptr == '^' && *(l_ptr + 1) == '(' && *(l_ptr + 2) == '-') {
            if (*(l_ptr - 1) == ')') {
              char const* l_ptr_c = l_ptr;

              --l_ptr_c;
              uint32_t parenthesis_counter = 0;

              while (*(l_ptr_c) != '(') {
                if (*l_ptr_c != ')' && *l_ptr_c != '(')
                  counter ++;
                else if (*l_ptr_c == ')')
                  parenthesis_counter++;

                l_ptr_c --;

                if ((is_operator(*l_ptr_c) || *l_ptr_c == '(') && is_operand(*(l_ptr_c + 1))) {
                  uint32_t tmp_c = 0;

                  while (is_operand(*(l_ptr_c + tmp_c + 2))) {
                    tmp_c ++;
                  }

                  counter -= tmp_c;
                } else if ((is_operator(*l_ptr_c) || *l_ptr_c == '(') && is_variable(*(l_ptr_c + 1))) {
                  uint32_t tmp_c = 0;

                  while (is_variable(*(l_ptr_c + tmp_c + 2)) || is_operand(*(l_ptr_c + tmp_c + 2))) {
                    tmp_c ++;
                  }

                  counter -= tmp_c;
                } else if (*l_ptr_c == '-' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                  counter += 1;
                else if (*l_ptr_c == '+' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                  counter += 1;
                else if (is_variable(*l_ptr_c) && *(l_ptr_c - 1) == '-' && *(l_ptr_c - 2) == '(')
                  counter -= 1;
                else if (is_variable(*l_ptr_c) && *(l_ptr_c - 1) == '+' && *(l_ptr_c - 2) == '(')
                  counter -= 1;
                else if (is_operand(*l_ptr_c) && *(l_ptr_c - 1) == '-' && *(l_ptr_c - 2) == '(')
                  counter -= 1;
                else if (is_operand(*l_ptr_c) && *(l_ptr_c - 1) == '+' && *(l_ptr_c - 2) == '(')
                  counter -= 1;

                if (*l_ptr_c == '(' && parenthesis_counter != 0) {
                  parenthesis_counter --;

                  if (parenthesis_counter == 0)
                    break;

                  l_ptr_c --;

                  if (*l_ptr_c == '^' && *(l_ptr_c + 1) == '(' && *(l_ptr_c + 2) == '-')
                    counter += 2; // one for '/' and one for '1'
                  else if (*l_ptr_c == '-' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                    counter += 1;
                  else if (*l_ptr_c == '+' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                    counter += 1;
                }

                if (parenthesis_counter == 0)
                  break;
              }
            } else
              counter = 1;

            neg_pow = true;
          } else if (*l_ptr == '^' && *(l_ptr + 1) == '-') {
            ERROR_MSG("Please put negative exponents in parentheses");
            std::exit(EXIT_FAILURE);
          }

          bool skip = false;

          if (neg_pow && *l_ptr == '^')
            op_stack.push('/');
          else if (!neg_pow && *l_ptr == '^' && pf.back().size() > 1 && pf.back()[0] == '-') {
            op_stack.push('!');
            pf.back().erase(pf.back().begin());
            skip = true;
          }

          if (!skip)
            op_stack.push(*l_ptr);
        }
      }
      // Push all open parenthesis to the stack
      else if (*l_ptr == '(')
        op_stack.push(*l_ptr);
      // When reaching a closing one, pop off operators from the stack until an opening one is found
      else if (*l_ptr == ')') {
        if (tmp.length() > 0) {
          if (neg_pow) {
            pf.insert(pf.end() - counter, "1");
            tmp.erase(0, 1);
            pf.emplace_back(tmp);
            tmp = "";
            counter = 0;
            neg_pow = false;
          } else {
            pf.emplace_back(tmp);
            tmp = "";
          }
        }

        tmp = "";

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

    if (tmp.length() > 0) {
      if (neg_pow) {
        pf.insert(pf.end() - counter, "1");
        tmp.erase(0, 1);
        pf.emplace_back(tmp);
        tmp = "";
        counter = 0;
        neg_pow = false;
      } else {
        pf.emplace_back(tmp);
      }
    }

    while (!op_stack.empty()) {
      pf.emplace_back(std::string(1, op_stack.top()));
      op_stack.pop();
    }

//      for (const auto & el : pf) {
//        std::cout << el << " ";
//      }
//
//      std::cout << "\n";
    pf.shrink_to_fit();
    functions.emplace_back(pf);
  }

  void ShuntingYardParser::parse_function(const std::string& fun, const std::vector<std::string>& vars, bool validate_fun) {
    std::string fun_ = fun;

    if (vars_map.empty()) {
      for (uint32_t i = 0; i != vars.size(); ++i) {
        vars_map.emplace(std::make_pair(vars[i], i));
      }
    }

    if (validate_fun)
      fun_ = validate(fun_, 0);

    parse(fun);

    functions.shrink_to_fit();
  }

  FFInt ShuntingYardParser::evaluate(const std::vector<std::string>& fun, const std::vector<FFInt>& values) const {
    FFInt res;

    std::stack<FFInt> nums;

    for (const auto& token : fun) {
      if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^" || token == "!") {
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

          nums.push(FFInt(mpz_class(tmp)));
        } else {
          if (token[0] == '-') {
            std::string tmp = token;
            tmp.erase(0, 1);

            if (vars_map.find(tmp) != vars_map.end())
              nums.push(-values[vars_map.at(tmp)]);
            else if (isdigit(tmp[0]))
              nums.push(-FFInt(std::stoull(tmp)));
            else
              throw_not_declared_var_err(tmp);
          } else if (token[0] == '+') {
            std::string tmp = token;
            tmp.erase(0, 1);

            if (vars_map.find(tmp) != vars_map.end())
              nums.push(values[vars_map.at(tmp)]);
            else  if (isdigit(tmp[0]))
              nums.push(FFInt(std::stoull(tmp)));
            else
              throw_not_declared_var_err(tmp);
          } else {
            if (vars_map.find(token) != vars_map.end())
              nums.push(values[vars_map.at(token)]);
            else if (isdigit(token[0]))
              nums.push(FFInt(std::stoull(token)));
            else
              throw_not_declared_var_err(token);
          }
        }
      }
    }

    if (nums.size())
      res = nums.top();
    else {
      ERROR_MSG("Error in functional evaluation! Please check your input.");
      std::exit(EXIT_FAILURE);
    }

    return res;
  }

  int ShuntingYardParser::get_weight(const char c) const {
    switch (c) {
      case '^':
      case '!':
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
    if (!is_operator(c) && c != '(' && c != ')' && chars.find(c) == chars.end() && c != ' ')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_operator(const char c) const {
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^')
      return true;

    return false;
  }

  bool ShuntingYardParser::is_variable(const char c) const {
    if (chars.find(c) != chars.end())
      return true;

    return false;
  }

  std::vector<std::vector<std::string>> ShuntingYardParser::get_rp_functions() const {
    return functions;
  }

  bool ShuntingYardParser::empty() const {
    return functions.empty();
  }

  void ShuntingYardParser::precompute_tokens() {
    precomp_tokens.clear();
    uint64_t size = functions.size();
    precomp_tokens = std::vector<std::vector<std::pair<uint8_t, FFInt>>> (size);

    for (uint64_t i = 0; i != size; ++i) {
      std::vector<std::string> tokens = functions[i];
      uint64_t t_size = tokens.size();
      precomp_tokens[i] = std::vector<std::pair<uint8_t, FFInt>> (t_size);

      uint32_t offset = 0;

      for (uint64_t j = 0; j != t_size; ++j) {
        std::string token = tokens[j + offset];

        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^" || token == "!") {
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
              // If the coefficient is a rational number, perform the division just once
              if (precomp_tokens[i][j - 1].first == operands::NUMBER && precomp_tokens[i][j - 2].first == operands::NUMBER) {
                FFInt quotient = precomp_tokens[i][j - 2].second / precomp_tokens[i][j - 1].second;
                j -= 2;
                t_size -= 2;
                offset += 2;
                precomp_tokens[i].pop_back();
                precomp_tokens[i].pop_back();
                precomp_tokens[i][j] = {operands::NUMBER, quotient};
              } else
                precomp_tokens[i][j] = {operands::OPERATOR, operators::DIV};

              break;
            }

            case '^': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::POW};
              break;
            }

            case '!': {
              precomp_tokens[i][j] = {operands::OPERATOR, operators::POW_NEG};
              break;
            }
          }
        } else {
          // check then if number has more than 18 digits
          if (token.length() > 18) {
            std::string tmp = token;

            if (token[0] == '+')
              tmp.erase(0, 1);

            precomp_tokens[i][j] = {operands::NUMBER, FFInt(mpz_class(tmp))};
          } else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {operands::NEG_VARIABLE, vars_map[tmp]};
              else if (std::isdigit(tmp[0]))
                precomp_tokens[i][j] = {operands::NUMBER, (-FFInt(std::stoull(tmp)))};
              else
                throw_not_declared_var_err(tmp);
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map[tmp]};
              else if (std::isdigit(tmp[0]))
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(tmp)))};
              else
                throw_not_declared_var_err(tmp);
            } else {
              if (vars_map.find(token) != vars_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map[token]};
              else if (std::isdigit(token[0]))
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(token)))};
              else
                throw_not_declared_var_err(token);
            }
          }
        }
      }
    }
  }

  std::string ShuntingYardParser::validate(const std::string& line, uint64_t exp_n) {
    size_t size = line.size() + 1;
    std::string r = line;
    std::stack<int> st;
    int i = 0;

    while (static_cast<size_t>(i) != size) {
      if (r[i] == '+' && r[i + 1] == '-')
        r[i] = '$';

      if (r[i] == '-' && r[i + 1] == '+')
        r[i + 1] = '$';

      if (r[i] == '(') {
        if (i != 0 && r[i - 1] == '(')
          st.push(-i);
        else
          st.push(i);

        i++;
      } else if (r[i] != ')' && r[i] != '(')
        i++;
      else if (r[i] == ')') {
        if (st.size() == 0) {
          ERROR_MSG("Mismatch of closing parentheses in expression " + std::to_string(exp_n) + ".");
          std::exit(EXIT_FAILURE);
        } else {
          int top = st.top();

          if (static_cast<uint32_t>(i) != size - 1 && r[i + 1] == ')' && top < 0) {
            r[-top] = '$';
            r[i] = '$';
          }

          st.pop();
          i++;
        }
      }
    }

    if (st.size() > 0) {
      ERROR_MSG("Mismatch of opening parentheses in expression " + std::to_string(exp_n) + ".");
      std::exit(EXIT_FAILURE);
    }

    std::string result = "";

    for (size_t tmp = 0; tmp != size; ++tmp) {
      if (r[tmp] == '$')
        continue;

      result += r[tmp];
    }

    result.shrink_to_fit();
    //std::cout << "validate\n" << line << "\n" << result << "\n";
    return result;
  }

  void ShuntingYardParser::throw_not_declared_var_err(const std::string& var) const {
    ERROR_MSG("Variable " + var + " not declared!");
    std::exit(EXIT_FAILURE);
  }
}
