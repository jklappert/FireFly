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

#include "firefly/HornerGenerator.hpp"
#include "firefly/Logger.hpp"

#include <map>

namespace firefly {

  std::string generate_horner_ff(const ff_map& monomials, const std::vector<std::string>& vars, uint32_t index) {
    if (!monomials.empty()) {
      std::map<uint32_t, ff_map, std::greater<uint32_t>> tmp_coefs {};

      if (monomials.begin() -> first.size() > 1) {
        for (const auto & mon : monomials) {
          uint32_t deg = mon.first[0];
          // Erase first entry to promote it to a monomial with n - 1 variables
          std::vector<uint32_t> degs = mon.first;
          degs.erase(degs.begin());

          if (tmp_coefs.find(deg) != tmp_coefs.end())
            tmp_coefs[deg].emplace(std::make_pair(degs, mon.second));
          else
            tmp_coefs.emplace(std::make_pair(deg, ff_map( {{degs, mon.second}})));
        }

        std::unordered_map<uint32_t, std::string> horner_coefs {};

        for (const auto & el : tmp_coefs) {
          horner_coefs.emplace(std::make_pair(el.first, generate_horner_ff(el.second, vars, index + 1)));
        }

        uint32_t max_deg = tmp_coefs.begin() -> first;
        std::string horner_coef = "";

        if (max_deg > 0) {
          for (uint32_t i = 0; i < max_deg - 1; ++i) {
            horner_coef += "(";
          }

          const std::string var = vars[index];

          horner_coef += var + "*(" + horner_coefs[max_deg] + ")";

          for (uint32_t i = max_deg - 1; i > 0; i--) {
            if (horner_coefs.find(i) != horner_coefs.end())
              horner_coef += "+" + horner_coefs[i];

            horner_coef += ")*" + var;
          }

          if (horner_coefs.find(0) != horner_coefs.end())
            horner_coef += "+" + horner_coefs[0];
        } else {
          if (horner_coefs.find(0) != horner_coefs.end())
            horner_coef += horner_coefs[0];
        }

        return horner_coef;
      } else if (monomials.begin() -> first.size() == 1) {
        uint32_t max_deg = 0;

        for (const auto & el : monomials) {
          uint32_t tmp_deg = el.first[0];

          if (tmp_deg > max_deg)
            max_deg = tmp_deg;
        }

        std::string horner_coef = "";

        if (max_deg > 0) {
          for (uint32_t i = 0; i < max_deg - 1; ++i) {
            horner_coef += "(";
          }

          const std::string var = vars[index];

          horner_coef += var;
          horner_coef += monomials.at( {max_deg}).n != 1 ? "*" + std::to_string(monomials.at( {max_deg}).n) : "";

          for (uint32_t i = max_deg - 1; i > 0; i--) {
            if (monomials.find( {i}) != monomials.end())
              horner_coef += "+" + std::to_string(monomials.at( {i}).n);

            horner_coef += ")*" + var;
          }

          if (monomials.find( {0}) != monomials.end())
            horner_coef += "+" + std::to_string(monomials.at( {0}).n);
        } else {
          if (monomials.find( {0}) != monomials.end())
            horner_coef += std::to_string(monomials.at( {0}).n);
        }

        return horner_coef;
      } else
        return std::to_string((monomials.begin() -> second).n);
    } else {
      WARNING_MSG("Provided an empty polynomial for Horner form. Will be interpreted as zero.");
      return "0";
    }
  }

  std::string generate_horner_rn(const rn_map& monomials, const std::vector<std::string>& vars, uint32_t index) {
    if (!monomials.empty()) {
      std::map<uint32_t, rn_map, std::greater<uint32_t>> tmp_coefs {};

      if (monomials.begin() -> first.size() > 1) {
        for (const auto & mon : monomials) {
          uint32_t deg = mon.first[0];
          // Erase first entry to promote it to a monomial with n - 1 variables
          std::vector<uint32_t> degs = mon.first;
          degs.erase(degs.begin());

          if (tmp_coefs.find(deg) != tmp_coefs.end())
            tmp_coefs[deg].emplace(std::make_pair(degs, mon.second));
          else
            tmp_coefs.emplace(std::make_pair(deg, rn_map( {{degs, mon.second}})));
        }

        std::unordered_map<uint32_t, std::string> horner_coefs {};

        for (const auto & el : tmp_coefs) {
          horner_coefs.emplace(std::make_pair(el.first, generate_horner_rn(el.second, vars, index + 1)));
        }

        uint32_t max_deg = tmp_coefs.begin() -> first;
        std::string horner_coef = "";

        if (max_deg > 0) {
          const std::string var = vars[index];

	  std::vector<uint32_t> non_zero_degrees;
	  non_zero_degrees.reserve(monomials.size());

          for (uint32_t i = max_deg; i > 0; i--) {
	    if (horner_coefs.find(i) != horner_coefs.end()) {
	      non_zero_degrees.emplace_back(i);
	      if (i != max_deg) {
		horner_coef += "(";
	      }
	    }
	  }

	  for (uint32_t i = 0; i != non_zero_degrees.size(); ++i) {
	    const uint32_t tmp_deg = non_zero_degrees[i];
	    if (i + 1 != non_zero_degrees.size()) {
	      const uint32_t h_deg = non_zero_degrees[i] - non_zero_degrees[i + 1];
	      if (tmp_deg == max_deg) {
		horner_coef += "(" + horner_coefs[tmp_deg] + ")*" + var;
		horner_coef += h_deg != 1 ? "^" + std::to_string(h_deg) : "";
	      } else {
		horner_coef += "+" + horner_coefs[tmp_deg] + ")*" + var;
		horner_coef += h_deg != 1 ? "^" + std::to_string(h_deg) : "";
	      }
	    } else {
	      const uint32_t h_deg = non_zero_degrees[i];
	      if (tmp_deg == max_deg) {
	        horner_coef += "(" + horner_coefs[tmp_deg] + ")*" + var;
		horner_coef += h_deg != 1 ? + "^" + std::to_string(h_deg) : "";
	      } else {
		horner_coef += "+" + horner_coefs[tmp_deg] + ")*" + var;
		horner_coef += h_deg != 1 ? + "^" + std::to_string(h_deg) : "";
	      }
	    }
	  }

          if (horner_coefs.find(0) != horner_coefs.end())
            horner_coef += "+" + horner_coefs[0];
        } else if (horner_coefs.find(0) != horner_coefs.end()) {
          horner_coef += horner_coefs[0];
        }

        return horner_coef;
      } else if (monomials.begin() -> first.size() == 1) {
        uint32_t max_deg = 0;

        for (const auto & el : monomials) {
          uint32_t tmp_deg = el.first[0];

          if (tmp_deg > max_deg)
            max_deg = tmp_deg;
        }

        std::string horner_coef = "";

        if (max_deg > 0) {
          const std::string var = vars[index];

	  std::vector<uint32_t> non_zero_degrees;
	  non_zero_degrees.reserve(monomials.size());

          for (uint32_t i = max_deg; i > 0; i--) {
	    if (monomials.find({i}) != monomials.end()) {
	      non_zero_degrees.emplace_back(i);
	      if (i != max_deg) {
		horner_coef += "(";
	      }
	    }
	  }

	  for (uint32_t i = 0; i != non_zero_degrees.size(); ++i) {
	    const uint32_t tmp_deg = non_zero_degrees[i];
	    if (i + 1 != non_zero_degrees.size()) {
	      const uint32_t h_deg = non_zero_degrees[i] - non_zero_degrees[i + 1];
	      if (tmp_deg == max_deg) {
		horner_coef += monomials.at({max_deg}).numerator != 1 || monomials.at({max_deg}).denominator != 1 ? monomials.at({tmp_deg}).string() + "*" + var :  var;
		horner_coef += h_deg != 1 ? "^" + std::to_string(h_deg) : "";
	      } else {
		horner_coef += "+" + monomials.at({tmp_deg}).string() + ")*" + var;
		horner_coef += h_deg != 1 ? "^" + std::to_string(h_deg) : "";
	      }
	    } else {
	      const uint32_t h_deg = non_zero_degrees[i];
	      if (tmp_deg == max_deg) {
	        horner_coef += monomials.at({max_deg}).numerator != 1 || monomials.at({max_deg}).denominator != 1 ? monomials.at({tmp_deg}).string() + "*" + var : var;
		horner_coef += h_deg != 1 ? + "^" + std::to_string(h_deg) : "";
	      } else {
		horner_coef += "+" + monomials.at({tmp_deg}).string() + ")*" + var;
		horner_coef += h_deg != 1 ? + "^" + std::to_string(h_deg) : "";
	      }
	    }
	  }

          if (monomials.find({0}) != monomials.end())
            horner_coef += "+" + monomials.at({0}).string();
        } else if (monomials.find({0}) != monomials.end()) {
          horner_coef += monomials.at({0}).string();//TODO do not add if coefficient is 1
        }

        return horner_coef;
      } else
        return (monomials.begin() -> second).string();
    } else {
      WARNING_MSG("Provided an empty polynomial to generate a Horner form. Will be interpreted as zero.");
      return "0";
    }
  }

  std::string generate_horner_mon(const std::vector<Monomial> monomials, const std::vector<std::string>& vars, uint32_t index) {
    if (!monomials.empty()) {
      rn_map map {};

      for (const auto & el : monomials) {
        map.emplace(el.powers, el.coef);
      }

      return generate_horner_rn(map, vars, index);
    } else {
      WARNING_MSG("Provided an empty polynomial for Horner form. Will be interpreted as zero.");
      return "0";
    }
  }
}
