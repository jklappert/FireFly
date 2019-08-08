#include "ShuntingYardParser.hpp"
#include "Logger.hpp"
#include <fstream>
#include <chrono>

namespace firefly {

  ShuntingYardParser::ShuntingYardParser() {}

  ShuntingYardParser::ShuntingYardParser(const std::string& file, const std::vector<std::string>& vars) {
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
    }

    uint32_t line_c = 1;

    while (std::getline(istream, line)) {
      line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

      if (line.length() > 0) {
        line = validate(line, line_c);
        parse(line);
      }

      line_c++;
    }

    functions.shrink_to_fit();
    precompute_tokens();

    istream.close();
    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Parsed " + std::to_string(functions.size()) + " functions in "
             + std::to_string(std::chrono::duration<double>(time1 - time0).count())
             + " s.");
  }

  void ShuntingYardParser::parse(const std::string& fun_) {
    std::string fun = fun_;

    // Check for global signs
    if (fun.size() > 2 && ((fun[0] == '+' || fun[0] == '-') && fun[1] == '(')) {
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
                  l_ptr_c --;

                  if (*l_ptr_c == '^' && *(l_ptr_c + 1) == '(' && *(l_ptr_c + 2) == '-')
                    counter += 2; // one for '/' and one for '1'
                  else if (*l_ptr_c == '-' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                    counter += 1;
                  else if (*l_ptr_c == '+' && *(l_ptr_c + 1) == '(' && *(l_ptr_c - 1) == '(')
                    counter += 1;

                  if (parenthesis_counter == 0)
                    break;
                }

                if (parenthesis_counter == 0)
                  break;
              }
            } else
              counter = 1;

            neg_pow = true;
          }

          if (neg_pow && *l_ptr == '^')
            op_stack.push('/');

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

//     for (const auto & el : pf) {
//       std::cout << el << " ";
//     }
//
//     std::cout << "\n";
    pf.shrink_to_fit();
    functions.emplace_back(pf);
  }

  void ShuntingYardParser::parse_function(const std::string& fun, const std::vector<std::string>& vars) {
    for (uint32_t i = 0; i < vars.size(); ++i) {
      vars_map.emplace(std::make_pair(vars[i], i));
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
          // check then if number has more than 18 digits
          if (token.length() > 18)
            nums.push(FFInt(mpz_class(token)));
          else {
            if (token[0] == '-') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                nums.push(-values[vars_map.at(tmp)]);
              else
                nums.push(-FFInt(std::stoull(tmp)));
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                nums.push(values[vars_map.at(tmp)]);
              else
                nums.push(FFInt(std::stoull(tmp)));
            } else {
              if (vars_map.find(token) != vars_map.end())
                nums.push(values[vars_map.at(token)]);
              else
                nums.push(FFInt(std::stoull(token)));
            }
          }
        }
      }

      if (nums.size())
        res.emplace_back(nums.top());
      else {
        ERROR_MSG("Error in functional evaluation! Please check your input.");
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
        ERROR_MSG("Error in functional evaluation! Please check your input.");
        std::exit(-1);
      }
    }

    return res;
  }

  std::vector<std::vector<FFInt>> ShuntingYardParser::evaluate_pre(const std::vector<std::vector<FFInt>>& values) const {
    size_t bunch_size = values.size();
    std::vector<std::vector<FFInt>> res(bunch_size);

    for (const auto & tokens : precomp_tokens) {
      std::vector<std::stack<FFInt>> nums(bunch_size);

      for (const auto & token : tokens) {
        switch (token.first) {
          case operands::OPERATOR : {
            // Pop two numbers
            for (size_t i = 0; i != bunch_size; ++i) {
              FFInt a = nums[i].top();
              nums[i].pop();
              FFInt b = nums[i].top();
              nums[i].pop();

              switch (token.second.n) {
                case operators::PLUS: {
                  nums[i].push(a + b);

                  break;
                }

                case operators::MINUS: {
                  nums[i].push(b - a);

                  break;
                }

                case operators::MULT: {
                  nums[i].push(b * a);

                  break;
                }

                case operators::DIV: {
                  nums[i].push(b / a);

                  break;
                }

                case operators::POW: {
                  nums[i].push(b.pow(a));

                  break;
                }
              }
            }

            break;
          }

          case operands::VARIABLE : {
            for (size_t i = 0; i != bunch_size; ++i) {
              nums[i].push(values[i][token.second.n]);
            }

            break;
          }

          case operands::NEG_VARIABLE : {
            for (size_t i = 0; i != bunch_size; ++i) {
              nums[i].push(-values[i][token.second.n]);
            }

            break;
          }

          case operands::NUMBER: {
            for (size_t i = 0; i != bunch_size; ++i) {
              nums[i].push(token.second);
            }
          }
        }
      }

      for (size_t i = 0; i != bunch_size; ++i) {
        if (nums[i].size())
          res[i].emplace_back(nums[i].top());
        else {
          ERROR_MSG("Error in functional evaluation! Check your input.");
          std::exit(-1);
        }
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

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {operands::NEG_VARIABLE, vars_map[tmp]};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (-FFInt(std::stoull(tmp)))};
            } else if (token[0] == '+') {
              std::string tmp = token;
              tmp.erase(0, 1);

              if (vars_map.find(tmp) != vars_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map[tmp]};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(tmp)))};
            } else {
              if (vars_map.find(token) != vars_map.end())
                precomp_tokens[i][j] = {operands::VARIABLE, vars_map[token]};
              else
                precomp_tokens[i][j] = {operands::NUMBER, (FFInt(std::stoull(token)))};
            }
          }
        }
      }
    }
  }

  std::string ShuntingYardParser::validate(const std::string& line, uint32_t exp_n) {
    size_t size = line.size() + 1;
    std::string r = line;
    std::stack<int> st;
    int i = 0;

    while (i < size) {
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
          ERROR_MSG("Mismatch of closing prenthesis in expression " + std::to_string(exp_n) + ".");
          std::exit(1);
        } else {
          int top = st.top();

          if (i != size - 1 && r[i + 1] == ')' && top < 0) {
            r[-top] = '$';
            r[i] = '$';
          }

          st.pop();
          i++;
        }
      }
    }

    if (st.size() > 0) {
      ERROR_MSG("Mismatch of opening prenthesis in expression " + std::to_string(exp_n) + ".");
      std::exit(1);
    }

    std::string result = "";

    for (i = 0; i < size; i++) {
      if (r[i] == '$')
        continue;

      result += r[i];
    }

    //std::cout << "validate\n" << line << "\n" << result << "\n";
    return result;
  }
}
