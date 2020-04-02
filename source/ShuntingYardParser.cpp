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
          validate(line, line_c);
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

        line_c++;
      }

      istream.close();
    } else {
      size_t curr_percentage = 1;

      for (size_t i = 0; i != funs.size(); ++i) {
        std::string fun = funs[i];
        validate(fun, i);
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

        if (i + 1 < funs.size() + 1 && i + 1 > curr_percentage*fun.size()/10) {
	  ++curr_percentage;
          std::cerr << "\033[1;34mFireFly info:\033[0m " << i + 1 << " / " << funs.size() << "\r";
	}
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

  void ShuntingYardParser::parse(std::string& fun) {
    //    std::string fun = fun_;

    // Check for global signs
    if (fun.size() > 2 && (fun[0] == '+' || fun[0] == '-') && fun[1] == '(') {
      if (fun[0] == '+')
        fun.erase(0, 1);
      else
        fun.insert(fun.begin(), '0');
    }

    std::string tmp = ""; // Used for numbers
    std::vector<std::string> tokens = {};
    std::stack<char> op_stack;
    size_t c_counter = 0;
    bool neg_exp = false;

    // Pick one character at a time until we reach the end of the line
    for (const char ex : fun) {
      // If operand, add it to postfix string
      // If operator pop operators off the stack until it is empty
      if (ex != '\0') {
        if (is_operand(ex)) {
          tmp.push_back(ex);
        } else if (is_variable(ex)) {
          tmp.push_back(ex);
        } else if (is_operator(ex)) {
          if (!neg_exp && tmp.length() > 0) {
            tokens.emplace_back(tmp);
            tmp = "";
          }

          // Check for cases like +(-(x+...))
          if (!op_stack.empty() && fun[c_counter - 1] == '(') {
            if (fun[c_counter + 1] == '(') {
              tokens.emplace_back("0");
              op_stack.push(ex);
            } else if (!neg_exp) {
              tmp.insert(tmp.begin(), ex);
            } else if (neg_exp) {
              neg_exp = false;
            }
          } else if (op_stack.empty() && tokens.empty()) {
            tmp.insert(tmp.begin(), ex);
	  } else {
            while (!op_stack.empty() && op_stack.top() != '(' && get_weight(op_stack.top()) >= get_weight(ex)) {
              tokens.emplace_back(std::string(1, op_stack.top()));
              op_stack.pop();
            }

            if (ex == '^' && fun[c_counter + 1] == '(' && fun[c_counter + 2] == '-') {
              neg_exp = true;
	      if (tokens.back().size() > 1 && tokens.back().front() == '-' && fun[c_counter - 1] != ')') {
		tokens.back().erase(tokens.back().begin());
		op_stack.push(';');
	      } else {
                op_stack.push('~');
	      }
            } else if (ex == '^' && tokens.back().size() > 1 && tokens.back().front() == '-' && fun[c_counter - 1] != ')') {
              op_stack.push('!');
              tokens.back().erase(tokens.back().begin());
            } else {
              op_stack.push(ex);
            }
          }
        }

        // Push all open parenthesis to the stack
        else if (ex == '(')
          op_stack.push(ex);
        // When reaching a closing one, pop off operators from the stack until an opening one is found
        else if (ex == ')') {
          if (tmp.length() > 0) {
            tokens.emplace_back(tmp);
            tmp = "";
          }

          tmp = "";

          while (!op_stack.empty()) {
            if (op_stack.top() == '(') {
              op_stack.pop();
              break;
            }

            tokens.emplace_back(std::string(1, op_stack.top()));
            op_stack.pop();
          }
        }

        ++c_counter;
      }
    }

    if (tmp.length() > 0) {
      tokens.emplace_back(tmp);
    }

    while (!op_stack.empty()) {
      tokens.emplace_back(std::string(1, op_stack.top()));
      op_stack.pop();
    }

    /*for (const auto& el : tokens) {
      std::cout << el << " ";
    }

    std::cout << "\n";*/

    tokens.shrink_to_fit();
    functions.emplace_back(tokens);
  }

  void ShuntingYardParser::parse_function(std::string& fun, const std::vector<std::string>& vars, bool validate_fun) {
    //    std::string fun_ = fun;

    if (vars_map.empty()) {
      for (uint32_t i = 0; i != vars.size(); ++i) {
        vars_map.emplace(std::make_pair(vars[i], i));
      }
    }

    if (validate_fun)
      validate(fun, 0);

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
      case '~':
      case ';':
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

        if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^" || token == "!" || token == "~" || token == ";") {
          switch (token[0]) {
            case '+': {
              precomp_tokens[i][j] = {tokens::PLUS, 0};
              break;
            }

            case '-': {
              precomp_tokens[i][j] = {tokens::MINUS, 0};
              break;
            }

            case '*': {
              precomp_tokens[i][j] = {tokens::MULT, 0};
              break;
            }

            case '/': {
              // If the coefficient is a rational number, perform the division just once
              if (precomp_tokens[i][j - 1].first == tokens::NUMBER && precomp_tokens[i][j - 2].first == tokens::NUMBER) {
                FFInt quotient = precomp_tokens[i][j - 2].second / precomp_tokens[i][j - 1].second;
                j -= 2;
                t_size -= 2;
                offset += 2;
                precomp_tokens[i].pop_back();
                precomp_tokens[i].pop_back();
                precomp_tokens[i][j] = {tokens::NUMBER, quotient};
              } else if (precomp_tokens[i][j - 1].first == tokens::NUMBER && (precomp_tokens[i][j - 2].first == tokens::VARIABLE || precomp_tokens[i][j - 2].first == tokens::NEG_VARIABLE)) {
                FFInt inverse = 1 / precomp_tokens[i][j - 1].second;
                precomp_tokens[i][j - 1].second = inverse;
                precomp_tokens[i][j] = {tokens::MULT, 0};
              } else
                precomp_tokens[i][j] = {tokens::DIV, 0};

              break;
            }

            case '^': {
              precomp_tokens[i][j] = {tokens::POW, 0};
              break;
            }

            case '~': {
              precomp_tokens[i][j] = {tokens::NEG_POW, 0};
              break;
            }

            case ';': {
              precomp_tokens[i][j] = {tokens::NEG_POW_NEG, 0};
              break;
            }

            case '!': {
              precomp_tokens[i][j] = {tokens::POW_NEG, 0};
              break;
            }
          }
        } else {
          // check then if number has more than 18 digits
          if (token.length() > 18) {
            std::string tmp = token;

            if (token[0] == '+')
              tmp.erase(0, 1);

            precomp_tokens[i][j] = {tokens::NUMBER, FFInt(mpz_class(tmp))};
          } else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {tokens::NEG_VARIABLE, vars_map[tmp]};
              else if (std::isdigit(tmp[0]))
                precomp_tokens[i][j] = {tokens::NUMBER, (-FFInt(std::stoull(tmp)))};
              else
                throw_not_declared_var_err(tmp);
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {tokens::VARIABLE, vars_map[tmp]};
              else if (std::isdigit(tmp[0]))
                precomp_tokens[i][j] = {tokens::NUMBER, (FFInt(std::stoull(tmp)))};
              else
                throw_not_declared_var_err(tmp);
            } else {
              if (vars_map.find(token) != vars_map.end())
                precomp_tokens[i][j] = {tokens::VARIABLE, vars_map[token]};
              else if (std::isdigit(token[0]))
                precomp_tokens[i][j] = {tokens::NUMBER, (FFInt(std::stoull(token)))};
              else
                throw_not_declared_var_err(token);
            }
          }
        }
      }
    }
  }

  void ShuntingYardParser::validate(std::string& r, uint64_t exp_n) {
    size_t size = r.size() + 1;
    //std::string r = line;
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

    r.erase(std::remove(r.begin(), r.end(), '$'), r.end());
    /*std::string result = "";

    for (size_t tmp = 0; tmp != size; ++tmp) {
      if (r[tmp] == '$')
        continue;

      result += r[tmp];
    }

    result.shrink_to_fit();
    r = std::move(result);*/
    //std::cout << "validate\n" << line << "\n" << result << "\n";
    //return result;
  }

  void ShuntingYardParser::throw_not_declared_var_err(const std::string& var) const {
    ERROR_MSG("Variable '" + var + "' not declared!");
    std::exit(EXIT_FAILURE);
  }
}
