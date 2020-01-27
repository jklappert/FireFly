//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
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

#include "Polynomial.hpp"
#include "HornerGenerator.hpp"
#include "Logger.hpp"
#include "utils.hpp"

namespace firefly {

  Polynomial::Polynomial(const rn_map& coef) {
    for (const auto & el : coef) {
      coefs.emplace_back(Monomial(el.first, el.second));
    }

    if (coefs.size() > 0)
      n = coefs[0].powers.size();
    else
      n = 0;
  }

  Polynomial::Polynomial(const Monomial& coef) {
    coefs.emplace_back(coef);
    n = coef.powers.size();
  }

  Polynomial::Polynomial() {}

  void Polynomial::sort() {
    std::sort(coefs.begin(), coefs.end(),
    [](const Monomial & a, const Monomial & b) {
      return a_grt_b(b.powers, a.powers);
    });
  }

  void Polynomial::clear() {
    coefs.clear();
  }

  std::string Polynomial::to_string(const std::vector<std::string>& vars) const {
    std::string str;

    if (coefs.empty())
      str += "0";
    else {
      if (vars.size() != n && var_pos == -1) {
        ERROR_MSG("Symbol size does not match to number of variables of the polynomial!");
        std::exit(EXIT_FAILURE);
      }

      std::vector<std::string> vars_c = vars;
      if (var_pos != -1) {
        vars_c[0] = vars[var_pos];
      }

      for (const auto & mono : coefs) {
        uint32_t deg = 0;
        for (const auto& el : mono.powers) {
          deg += el;
        }

        if (deg != 0 && (mono.coef.numerator != 1 || mono.coef.denominator != 1))
          str += mono.coef.string() + "*";
        else if (deg == 0)
          str += mono.coef.string() + "*";

        for (uint32_t i = 0; i != mono.powers.size(); ++i) {
          if (mono.powers[i] > 1) {
            str += vars_c[i] + "^" + std::to_string(mono.powers[i]) + "*";
          } else if (mono.powers[i] == 1) {
            str += vars_c[i] + "*";
          }
        }

        str.erase(--str.end());

        str += "+";
      }

      str.erase(--str.end());
    }

    return str;
  }

  std::ostream& operator<<(std::ostream& out, const Polynomial& pol) {
    bool first = true;

    for (const auto & mono : pol.coefs) {
      if (first) {
        out <<  mono.coef << "*x^(";

        for (const auto i : mono.powers) {
          out << i << ",";
        }

        out << "\b)";
        first = false;
      } else {
        out << " + " << mono.coef << "*x^(";

        for (const auto i : mono.powers) {
          out << i << ",";
        }

        out << "\b)";
      }

    }

    if (pol.coefs.empty())
      out << "0";

    out << "\n";
    return out;
  }

  PolynomialFF Polynomial::convert_to_PolynomialFF() {
    ff_map coefs_ff;
    uint32_t n = coefs[0].powers.size();

    for (auto & coef : coefs) {
      FFInt coef_ff = FFInt(coef.coef.numerator) / FFInt(coef.coef.denominator);

      if (coef_ff.n > 0)
        coefs_ff.emplace(std::make_pair(coef.powers, coef_ff));
    }

    return PolynomialFF(n, coefs_ff);
  }

  Polynomial& Polynomial::operator*=(const RationalNumber& rn) {
    for (auto & mon : coefs) {
      mon.coef *= rn;
    }

    return *this;
  }

  std::string Polynomial::generate_horner(std::vector<std::string> vars) const {
    return generate_horner_mon(coefs, vars);
  }

  void Polynomial::set_var_pos(int var_pos_) {
    var_pos = var_pos_;
  }
}
