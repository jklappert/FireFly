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

#include "firefly/gzstream.hpp"
#include "firefly/ParserUtils.hpp"
#include "firefly/tinydir.h"

namespace firefly {
  std::vector<uint32_t> parse_vector_32(std::string& line, int number_of_parameters) {
    size_t pos = 0;
    int i = 0;
    std::string delimiter = " ";
    std::vector<uint32_t> tmp {};

    if (number_of_parameters > 0)
      tmp.reserve(number_of_parameters);

    while ((pos = line.find(delimiter)) != std::string::npos) {
      tmp.emplace_back(std::stoi(line.substr(0, pos)));
      line.erase(0, pos + 1);
      ++i;

      if (i == number_of_parameters) break;
    }

    return tmp;
  }

  std::vector<FFInt> parse_vector_FFInt(std::string& line, int number_of_parameters) {
    size_t pos = 0;
    int i = 0;
    std::string delimiter = " ";
    std::vector<FFInt> tmp {};

    if (number_of_parameters > 0)
      tmp.reserve(number_of_parameters);

    while ((pos = line.find(delimiter)) != std::string::npos) {
      tmp.emplace_back(std::stoul(line.substr(0, pos)));
      line.erase(0, pos + 1);
      ++i;

      if (i == number_of_parameters) break;
    }

    return tmp;
  }

  RationalNumber parse_rational_number(const std::string& line) {
    size_t pos = line.find(" ");
    return RationalNumber(mpz_class(line.substr(0, pos)), mpz_class(line.substr(pos + 1)));
  }

  uint32_t parse_prime_number(const std::string& file_name) {
    std::string reverse_file_name = file_name;
    std::reverse(reverse_file_name.begin(), reverse_file_name.end());
    reverse_file_name.erase(0, 3);
    size_t pos = reverse_file_name.find("_");
    std::string prime = reverse_file_name.substr(0, pos);
    std::reverse(prime.begin(), prime.end());
    return std::stoi(prime);
  }

  std::unordered_map<uint32_t, std::list<RationalFunction>> parse_factors_rf() {
    std::unordered_map<uint32_t, std::list<RationalFunction>> factors_rf {};
    std::string line;
    tinydir_dir fac_dir;
    tinydir_open_sorted(&fac_dir, "ff_save/factors_rf");

    std::vector<std::string> fac_files;

    for (size_t i = 0; i != fac_dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&fac_dir, &file, i);

      if (!file.is_dir) {
        fac_files.emplace_back(file.name);
      }
    }

    tinydir_close(&fac_dir);

    for (const auto & file : fac_files) {
      std::string fac_number = "";
      for (const auto & character : file) {
        if (character != '.') {
          fac_number += character;
        } else {
          break;
        }
      }

      rn_map tmp_numerator;
      rn_map tmp_denominator;

      igzstream fac_file;
      std::string fac_path = "ff_save/factors_rf/" + file;
      fac_file.open(fac_path.c_str());
      bool is_num = false;
      bool is_var = true;
      int var_pos = -1;

      while (std::getline(fac_file, line)) {
        if (line == "var") {
          is_var = true;
          is_num = false;

          if (!tmp_numerator.empty() || !tmp_denominator.empty()) {
            Polynomial tmp_num(tmp_numerator);
            Polynomial tmp_den(tmp_denominator);
            tmp_num.set_var_pos(var_pos);
            tmp_den.set_var_pos(var_pos);
            factors_rf[std::stoi(fac_number)].emplace_back(RationalFunction(tmp_num, tmp_den));
          }

          tmp_numerator.clear();
          tmp_denominator.clear();
        } else if (line == "numerator") {
          is_num = true;
        } else if (line == "denominator") {
          is_num = false;
        } else if(is_var) {
          is_var = false;
          var_pos = std::stoi(line);
        } else if(is_num) {
          std::vector<uint32_t> tmp_vec = parse_vector_32(line, 1);
          tmp_numerator.emplace(std::make_pair(tmp_vec, parse_rational_number(line)));
        } else {
          std::vector<uint32_t> tmp_vec = parse_vector_32(line, 1);
          tmp_denominator.emplace(std::make_pair(tmp_vec, parse_rational_number(line)));
        }
      }

      if (!tmp_numerator.empty() || !tmp_denominator.empty()) {
        Polynomial tmp_num(tmp_numerator);
        Polynomial tmp_den(tmp_denominator);
        tmp_num.set_var_pos(var_pos);
        tmp_den.set_var_pos(var_pos);
        factors_rf[std::stoi(fac_number)].emplace_back(RationalFunction(tmp_num, tmp_den));
      }

      fac_file.close();
    }

    return factors_rf;
  }

  std::unordered_map<uint32_t, ShuntingYardParser> parse_factors(uint32_t n) {
    std::unordered_map<uint32_t, ShuntingYardParser> parsed_factors {};
    std::string line;
    tinydir_dir fac_dir;
    tinydir_open_sorted(&fac_dir, "ff_save/factors");

    std::vector<std::string> fac_files;

    for (size_t i = 0; i != fac_dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&fac_dir, &file, i);

      if (!file.is_dir) {
        fac_files.emplace_back(file.name);
      }
    }

    tinydir_close(&fac_dir);

    std::vector<std::string> fac_vars (n);
    for (size_t i = 0; i != n; ++i) {
      fac_vars[i] = "x" + std::to_string(i + 1);
    }

    for (const auto & file : fac_files) {
      std::string fac_number = "";
      for (const auto & character : file) {
        if (character != '.') {
          fac_number += character;
        } else {
          break;
        }
      }

      igzstream fac_file;
      std::string fac_path = "ff_save/factors/" + file;
      fac_file.open(fac_path.c_str());
      std::getline(fac_file, line);
      line.pop_back();
      ShuntingYardParser parser = ShuntingYardParser();
      parser.parse_function(line, fac_vars);
      parser.precompute_tokens();
      parsed_factors.emplace(std::stoi(fac_number), parser);
      fac_file.close();
    }

    return parsed_factors;
  }
}
