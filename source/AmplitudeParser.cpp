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

#include "AmplitudeParser.hpp"

#include <chrono>
#include <fstream>

namespace firefly {

  AmplitudeParser::AmplitudeParser() {}

  AmplitudeParser::AmplitudeParser(const std::vector<std::string>& vars_, const std::vector<std::string>& integral_families_) : vars(vars_), integral_families(integral_families_) {}

  void AmplitudeParser::parse_amplitude_file(const std::string& amplitude_file) {
    logger.open("ff_insert.log", std::ios_base::app);
    INFO_MSG("Parsing expression of " + amplitude_file);
    logger << "Parsing expression of " << amplitude_file << "\n";
    logger.close();
    parse_file(amplitude_file);
  }

  void AmplitudeParser::parse_file(const std::string& amplitude_file) {
    // Check if file exists
    std::ifstream infile_test(amplitude_file);
    std::ifstream infile;

    if (!infile_test.good()) {
      ERROR_MSG("File '" + amplitude_file + "' does not exist!");
      std::exit(EXIT_FAILURE);
    }

    infile.open(amplitude_file);

    std::string line;
    size_t counter = 0;

    while (std::getline(infile, line, ';')) {
      if (counter == 1)
        break;

      parse_amplitude_string(line);
      ++counter;
    }
  }

  void AmplitudeParser::parse_amplitude_string(const std::string& amplitude_string) {
    auto time0 = std::chrono::high_resolution_clock::now();
    auto res = parse_string(amplitude_string);
    functions.reserve(functions.size() + res.size());

    distinct_integral_counter = 0;
    parser_counter = 0;

    for (const auto& int_coef_pair : res) {
      if (integrals.find(int_coef_pair.first) != integrals.end()) {
        functions.emplace_back(int_coef_pair.second);
        amplitude_mapping[integrals[int_coef_pair.first]].emplace_back(std::make_pair(parser_counter, coef_type::PREFACTOR));
        ++parser_counter;
      } else {
        integrals.emplace(std::make_pair(int_coef_pair.first, distinct_integral_counter));
        functions.emplace_back(int_coef_pair.second);
        amplitude_mapping.emplace(std::make_pair(distinct_integral_counter, std::vector<std::pair<size_t, size_t>> (1, std::make_pair(parser_counter, coef_type::PREFACTOR))));
        ++parser_counter;
        ++distinct_integral_counter;
      }
    }

    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Parsed expression in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s");
    INFO_MSG("Found " + std::to_string(distinct_integral_counter) + " distinct functions\n");
    logger.open("ff_insert.log", std::ios_base::app);
    logger << "Parsed expression in " << std::to_string(std::chrono::duration<double>(time1 - time0).count()) << " s\n";
    logger << "Found " << std::to_string(distinct_integral_counter) << " distinct functions\n\n";
    logger.close();
  }

  std::vector<std::pair<std::string, std::string>> AmplitudeParser::parse_string(const std::string& amplitude) {
    std::vector<std::pair<std::string, std::string>> parsed_integrals {};

    std::string amplitude_c = amplitude;
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), ' '), amplitude_c.end());
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), '\t'), amplitude_c.end());
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), '\n'), amplitude_c.end());
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), '{'), amplitude_c.end());
    amplitude_c.erase(std::remove(amplitude_c.begin(), amplitude_c.end(), '}'), amplitude_c.end());
    std::size_t amplitude_size = amplitude_c.size();

    if (amplitude_c.length() > 0) {
      std::size_t old_pos = 0;
      std::size_t found = amplitude_c.find(']');

      while (found != std::string::npos) {
        std::string seed = amplitude_c.substr(old_pos, found + 1 - old_pos);
        std::size_t fo = seed.find('[');
        std::string int_fam = seed.substr(0, fo);

        bool found_integral = false;

        for (const auto& fam : integral_families) {
          size_t fo2 = int_fam.find(fam);
          std::string coefficient;

          if (fo2 != std::string::npos) {
            found_integral = true;
            std::string prefac = "";

            std::string integral = int_fam.substr(fo2) + seed.substr(fo);

            if (fo2 != 0) {
              prefac = int_fam.substr(fo2 - 1, 1);
            }

            if (prefac == "*") {
              coefficient = int_fam.substr(0, int_fam.size() - fam.size() - 1);
            } else if (found + 1 < amplitude_size - 1 && amplitude_c.substr(found + 1, 1) == "*") {
              std::size_t fo_2 = amplitude_c.find('[', found + 1);
	      if (fo_2 != std::string::npos) {
                std::string seed_2 = amplitude_c.substr(found + 1, fo_2-(found+1));
		for (const auto& tmp_fam : integral_families) {
		  std::size_t tmp_found = seed_2.rfind(tmp_fam);
		  if (tmp_found != std::string::npos) {
		    if (seed_2[0] == '*')
		      coefficient = seed_2.substr(1,tmp_found - 2);
		    else
		      coefficient = seed_2.substr(0,tmp_found - 2);

		    found = fo_2 - tmp_fam.size() - 2;
		  }
		}
              } else {
                coefficient = amplitude_c.substr(found + 2);
	      }
            } else
              coefficient = "1";

            if (prefac == "-")
              parsed_integrals.emplace_back(std::make_pair(integral, prefac + "(" + coefficient + ")"));
            else
              parsed_integrals.emplace_back(std::make_pair(integral, coefficient));
          }

          if (found_integral)
            break;
        }

        if (!found_integral) {
          ERROR_MSG("Unknown function: " + int_fam);
          logger.open("ff_insert.log", std::ios_base::app);
          logger << "Unknown function: " << int_fam << "\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }

        old_pos = found + 1;
        found = amplitude_c.find(']', found + 1);
      }
    } else {
      ERROR_MSG("Expression has no content");
      logger.open("ff_insert.log", std::ios_base::app);
      logger << "Expression has no content\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    parsed_integrals.shrink_to_fit();

    /*for (const auto& el : parsed_integrals) {
      std::cout << "Integral: " << el.first << "\n";
      std::cout << "Coefficient: " << el.second << "\n";
      }*/

    return parsed_integrals;
  }

  void AmplitudeParser::parse_ibp_table_file(const std::string& ibp_table) {
    // Check if file exists
    std::ifstream infile_test(ibp_table);
    std::ifstream infile;

    if (!infile_test.good()) {
      ERROR_MSG("File '" + ibp_table + "' does not exist!");
      logger.open("ff_insert.log", std::ios_base::app);
      logger << "File '" << ibp_table << "' does not exist!\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    INFO_MSG("Parsing replacement table of " + ibp_table);
    logger.open("ff_insert.log", std::ios_base::app);
    logger << "Parsing replacement table of " << ibp_table << "\n";
    logger.close();

    infile.open(ibp_table);

    std::string line;

    size_t counter = 0;

    while (std::getline(infile, line, '}')) {
      if (counter == 1)
        break;

      parse_ibp_table_string(line);
      ++counter;
    }
  }

  void AmplitudeParser::parse_ibp_table_string(const std::string& ibp_table) {
    auto time0 = std::chrono::high_resolution_clock::now();

    std::string ibp_table_c = ibp_table;
    ibp_table_c.erase(std::remove(ibp_table_c.begin(), ibp_table_c.end(), ' '), ibp_table_c.end());
    ibp_table_c.erase(std::remove(ibp_table_c.begin(), ibp_table_c.end(), '\t'), ibp_table_c.end());
    ibp_table_c.erase(std::remove(ibp_table_c.begin(), ibp_table_c.end(), '\n'), ibp_table_c.end());
    ibp_table_c.erase(std::remove(ibp_table_c.begin(), ibp_table_c.end(), '{'), ibp_table_c.end());
    ibp_table_c.erase(std::remove(ibp_table_c.begin(), ibp_table_c.end(), '}'), ibp_table_c.end());

    size_t pos = ibp_table_c.find("->");
    size_t old_pos = 0;
    size_t required_repl_counter = 0;

    while (pos != std::string::npos) {
      std::string lhs = ibp_table_c.substr(old_pos, pos - old_pos);
      std::string rhs;
      old_pos = pos + 2;
      size_t tmp_pos = ibp_table_c.find("->", old_pos);

      // reached end of table
      if (tmp_pos == std::string::npos) {
        rhs = ibp_table_c.substr(old_pos);
        auto repl_lhs = parse_string(lhs);

        if (integrals.find(repl_lhs[0].first) != integrals.end()) {
          if (repl_integrals.find(repl_lhs[0].first) != repl_integrals.end()) {
            ERROR_MSG("Multiple replacement rules for function: " + repl_lhs[0].first);
            logger.open("ff_insert.log", std::ios_base::app);
            logger << "Multiple replacement rules for function: " << repl_lhs[0].first << "\n";
            logger.close();
            std::exit(EXIT_FAILURE);
          }

          repl_integrals.emplace(repl_lhs[0].first);

          ++required_repl_counter;
          auto repl_rhs = parse_string(rhs);

          for (const auto & mi_coef : repl_rhs) {
            if (masters.find(mi_coef.first) == masters.end()) {
              masters.emplace(std::make_pair(mi_coef.first, distinct_master_counter));
              masters_inv.emplace(std::make_pair(distinct_master_counter - coef_type::COEF_TYPE_SIZE, mi_coef.first));
              ++distinct_master_counter;
            }

            amplitude_mapping[integrals[repl_lhs[0].first]].emplace_back(std::make_pair(parser_counter, masters[mi_coef.first]));
            functions.emplace_back(mi_coef.second);

            ++parser_counter;
          }
        }

        break;
      } else {
        tmp_pos = ibp_table_c.rfind("[", tmp_pos);
        tmp_pos = ibp_table_c.rfind(",", tmp_pos);
        rhs = ibp_table_c.substr(old_pos, tmp_pos - old_pos);
        old_pos = tmp_pos;
        pos = ibp_table_c.find("->", old_pos);
      }

      auto repl_lhs = parse_string(lhs);

      if (integrals.find(repl_lhs[0].first) != integrals.end()) {
        if (repl_integrals.find(repl_lhs[0].first) != repl_integrals.end()) {
          ERROR_MSG("Multiple replacement rules for function: " + repl_lhs[0].first);
          logger.open("ff_insert.log", std::ios_base::app);
          logger << "Multiple replacement rules for function: " << repl_lhs[0].first << "\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }

        repl_integrals.emplace(repl_lhs[0].first);

        ++required_repl_counter;
        auto repl_rhs = parse_string(rhs);

        for (const auto & mi_coef : repl_rhs) {
          if (masters.find(mi_coef.first) == masters.end()) {
            masters.emplace(std::make_pair(mi_coef.first, distinct_master_counter));
            masters_inv.emplace(std::make_pair(distinct_master_counter - coef_type::COEF_TYPE_SIZE, mi_coef.first));
            ++distinct_master_counter;
          }

          amplitude_mapping[integrals[repl_lhs[0].first]].emplace_back(std::make_pair(parser_counter, masters[mi_coef.first]));
          functions.emplace_back(mi_coef.second);

          ++parser_counter;
        }
      }
    }

    functions.shrink_to_fit();
    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Parsed replacement table in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s");
    INFO_MSG("Found " + std::to_string(required_repl_counter) + " required replacement rules");
    INFO_MSG("Found " + std::to_string(distinct_master_counter - coef_type::COEF_TYPE_SIZE) + " basis function(s) in total\n");
    logger.open("ff_insert.log", std::ios_base::app);
    logger << "Parsed replacement table in " << std::to_string(std::chrono::duration<double>(time1 - time0).count()) << " s\n";
    logger << "Found " << std::to_string(required_repl_counter) << " required replacement rules\n";
    logger << "Found " << std::to_string(distinct_master_counter - coef_type::COEF_TYPE_SIZE) << " basis function(s) in total\n\n";
    logger.close();
  }

  size_t AmplitudeParser::check_for_unreplaced_masters() {
    bool found_unreplaced_integral = false;

    for (const auto & integral : integrals) {
      if (repl_integrals.find(integral.first) == repl_integrals.end()) {
        found_unreplaced_integral = true;
        INFO_MSG("No replacements found for function: " + integral.first);
        logger.open("ff_insert.log", std::ios_base::app);
        logger << "No replacements found for function: " << integral.first << "\n";
        logger.close();

        if (masters.find(integral.first) == masters.end()) {
          masters.emplace(std::make_pair(integral.first, distinct_master_counter));
          masters_inv.emplace(std::make_pair(distinct_master_counter - coef_type::COEF_TYPE_SIZE, integral.first));
          amplitude_mapping[integral.second].emplace_back(std::make_pair(parser_counter, distinct_master_counter));
          ++distinct_master_counter;
        } else
          amplitude_mapping[integral.second].emplace_back(std::make_pair(parser_counter, masters[integral.first]));

        functions.emplace_back("1");
        ++parser_counter;
      }
    }

    if (found_unreplaced_integral) {
      INFO_MSG("Found " + std::to_string(distinct_master_counter - coef_type::COEF_TYPE_SIZE) + " basis function(s) in total\n");
      logger.open("ff_insert.log", std::ios_base::app);
      logger << "Found " << std::to_string(distinct_master_counter - coef_type::COEF_TYPE_SIZE) << " basis function(s) in total\n\n";
      logger.close();
    }

    return distinct_master_counter - coef_type::COEF_TYPE_SIZE;
  }

  std::string AmplitudeParser::get_unsimplified_coef(size_t master) const {
    std::string tmp_fun = "";

    for (size_t i = 0; i != distinct_integral_counter; ++i) {
      std::string tmp_coef = "+(";
      bool got_master = false;

      for (const auto& mi_map : amplitude_mapping.at(i)) {
        bool pref_done = false;

        if (mi_map.second == coef_type::PREFACTOR || mi_map.second == coef_type::REPEATED_REP)
          tmp_coef += "+(" + functions[mi_map.first] + ")";
        else {
          if (mi_map.second - coef_type::COEF_TYPE_SIZE == master) {
            if (!pref_done) {
              pref_done = true;
              tmp_coef += ")*(";
            }

            tmp_coef += "+(" + functions[mi_map.first] + ")";
            got_master = true;
          }
        }
      }

      tmp_coef += ")";

      if (got_master)
        tmp_fun += tmp_coef;
    }

    return tmp_fun;
  }

  FFAmplitudeBlackBox AmplitudeParser::build_black_box(size_t master) const {
    std::string tmp_fun = "";

    for (size_t i = 0; i != distinct_integral_counter; ++i) {
      std::string tmp_coef = "+(";
      bool got_master = false;

      for (const auto& mi_map : amplitude_mapping.at(i)) {
        bool pref_done = false;

        if (mi_map.second == coef_type::PREFACTOR || mi_map.second == coef_type::REPEATED_REP)
          tmp_coef += "+(" + functions[mi_map.first] + ")";
        else {
          if (mi_map.second - coef_type::COEF_TYPE_SIZE == master) {
            if (!pref_done) {
              pref_done = true;
              tmp_coef += ")*(";
            }

            tmp_coef += "+(" + functions[mi_map.first] + ")";
            got_master = true;
          }
        }
      }

      tmp_coef += ")";

      if (got_master)
        tmp_fun += tmp_coef;
    }

    return FFAmplitudeBlackBox(vars, {tmp_fun});
  }

  std::string AmplitudeParser::get_master(size_t i) {
    return masters_inv.at(i);
  }
}
