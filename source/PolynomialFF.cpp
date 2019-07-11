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

#include "PolynomialFF.hpp"
#include "Logger.hpp"
//#include <map>
//#include <chrono>

namespace firefly {

  PolynomialFF::PolynomialFF(uint32_t n_, ff_map coefs_) : n(n_), coefs(coefs_) {}

  PolynomialFF::PolynomialFF() {}

  FFInt PolynomialFF::calc(const std::vector<FFInt>& x) const {
    FFInt res(0);

    for (const auto & term : coefs) {
      FFInt product(1);

      for (uint32_t i = 0; i < n; ++i) {
        product *= x[i].pow(term.first[i]);
      }

      res += term.second * product;
    }

    return res;
  }

  FFInt PolynomialFF::calc_n_m_1(const std::vector<FFInt>& x) const {
    FFInt res(0);
    uint32_t n_m_1 = n - 1;

    for (const auto & term : coefs) {
      FFInt product(1);

      for (uint32_t i = 0; i < n_m_1; ++i) {
        product *= x[i].pow(term.first[i + 1]);
      }

      res += term.second * product;
    }

    return res;
  }

  std::unordered_map<uint32_t, FFInt> PolynomialFF::calc_n_m_1_map(const std::vector<FFInt>& x) const {
    std::unordered_map<uint32_t, FFInt> eval_map {};
    uint32_t n_m_1 = n - 1;

    for (const auto & term : coefs) {
      FFInt product(1);
      uint32_t deg = 0;

      for (uint32_t i = 0; i < n_m_1; ++i) {
        deg += term.first[i + 1];
        product *= x[i].pow(term.first[i + 1]);
      }

      deg += term.first[0];

      if (eval_map.find(deg) != eval_map.end())
        eval_map[deg] += term.second * product;
      else
        eval_map.emplace(std::make_pair(deg, term.second * product));
    }

    return eval_map;
  }

  PolynomialFF operator+(const PolynomialFF& a, const PolynomialFF& b) {
    ff_map new_coefs {};

    if (a.coefs.size() < b.coefs.size()) {
      new_coefs = b.coefs;

      for (const auto & el : a.coefs) {
        auto got = new_coefs.find(el.first);

        if (got == new_coefs.end()) {
          new_coefs.emplace(el);
        } else {
          got -> second += el.second;

          if (got -> second == 0)
            new_coefs.erase(el.first);
        }
      }
    } else {
      new_coefs = a.coefs;

      for (const auto & el : b.coefs) {
        auto got = new_coefs.find(el.first);

        if (got == new_coefs.end()) {
          new_coefs.emplace(el);
        } else {
          got -> second += el.second;

          if (got -> second == 0)
            new_coefs.erase(el.first);
        }
      }
    }

    return PolynomialFF(a.n, new_coefs);
  }

  PolynomialFF operator-(const PolynomialFF& a, const PolynomialFF& b) {
    ff_map new_coefs = a.coefs;

    for (const auto & el : b.coefs) {
      auto got = new_coefs.find(el.first);

      if (got == new_coefs.end()) {
        FFInt num = FFInt(0) - el.second;
        new_coefs.emplace(std::make_pair(el.first, num));
      } else {
        got -> second -= el.second;

        if (got -> second == 0)
          new_coefs.erase(el.first);
      }
    }

    return PolynomialFF(a.n, new_coefs);
  }

  PolynomialFF& PolynomialFF::operator-=(const PolynomialFF& b) {
    for (const auto & coef_b : b.coefs) {
      if (coef_b.second != 0) {
        if (coefs.find(coef_b.first) == coefs.end())
          coefs[coef_b.first] = -coef_b.second;
        else {
          coefs[coef_b.first] -= coef_b.second;

          if (coefs[coef_b.first] == 0)
            coefs.erase(coef_b.first);
        }
      }
    }

    return *this;
  }

  PolynomialFF& PolynomialFF::operator+=(const PolynomialFF& b) {
    for (const auto & coef_b : b.coefs) {
      if (coef_b.second != 0) {
        if (coefs.find(coef_b.first) == coefs.end())
          coefs[coef_b.first] = coef_b.second;
        else {
          coefs[coef_b.first] += coef_b.second;

          if (coefs[coef_b.first] == 0)
            coefs.erase(coef_b.first);
        }
      }
    }

    return *this;
  }

  PolynomialFF operator*(const PolynomialFF& a, const FFInt& ffint) {
    ff_map new_coefs = a.coefs;

    for (auto & el : new_coefs) {
      el.second *= ffint;
    }

    return PolynomialFF(a.n, new_coefs);
  }


  PolynomialFF operator/(const PolynomialFF& a, const FFInt& ffint) {
    ff_map new_coefs = a.coefs;
    FFInt inv = 1 / ffint;

    for (auto & el : new_coefs) {
      el.second *= inv;
    }

    return PolynomialFF(a.n, new_coefs);
  }

  PolynomialFF& PolynomialFF::operator*=(const FFInt& ffint) {
    for (auto & el : coefs) {
      el.second *= ffint;
    }

    return *this;
  }

  PolynomialFF& PolynomialFF::operator/=(const FFInt& ffint) {
    FFInt inv = 1 / ffint;

    for (auto & el : coefs) {
      el.second *= inv;
    }

    return *this;
  }

  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a) {
    for (auto & coef_ : a.coefs) {
      out << " + " << coef_.second.n << "*x^(";

      for (const auto i : coef_.first) {
        out << i << ",";
      }

      out << "\b)";
    }

    out << "\n";
    return out;
  }

  PolynomialFF PolynomialFF::mul(const uint32_t zi) const {
    ff_map new_coefs {};
    new_coefs.reserve(coefs.size());

    for (const auto & coef_ : coefs) {
      std::vector<uint32_t> new_element = coef_.first;
      ++new_element[zi - 1];
      new_coefs.emplace(std::make_pair(new_element, coef_.second));
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::homogenize(uint32_t degree) {
    ff_map tmp_coef = coefs;
    coefs.clear();

    for (const auto & mon : tmp_coef) {
      uint32_t old_degree = 0;

      std::vector<uint32_t> old_powers = mon.first;
      std::vector<uint32_t> new_powers = old_powers;

      for (const auto & power : old_powers) old_degree += power;

      new_powers.emplace(new_powers.begin(), degree - old_degree);
      coefs[new_powers] = tmp_coef[old_powers];
    }

    return *this;
  }

  bool PolynomialFF::zero() const {
    if (coefs.empty())
      return true;
    else if (coefs.size() == 1 && coefs.begin()->second == 0)
      return true;

    return false;
  }

  std::vector<uint32_t> PolynomialFF::max_deg() {
    if (max_degree.empty()) {
      int tmp_max;
      int tmp_min;

      for (const auto c : coefs) {
        int tmp_deg = 0;

        for (const auto i : c.first) {
          tmp_deg += i;
        }

        if (max_degree.empty()) {
          tmp_max = tmp_deg;
          tmp_min = tmp_deg;
        }

        tmp_max = std::max(tmp_deg, tmp_max);
        tmp_min = std::min(tmp_deg, tmp_min);

        if (tmp_max == tmp_deg) max_degree = c.first;

        if (tmp_min == tmp_deg) min_degree = c.first;
      }
    }

    return max_degree;
  }

  std::vector<uint32_t> PolynomialFF::min_deg() {
    max_deg();
    return min_degree;
  }

  PolynomialFF operator*(const PolynomialFF& a, const PolynomialFF& b) {
    ff_map new_monomials;
    new_monomials.reserve(a.coefs.size()*b.coefs.size() + 1);

    for (const auto & coef_a : a.coefs) {
      for (const auto & coef_b : b.coefs) {

        FFInt new_coef = coef_a.second * coef_b.second;

        if (new_coef != 0) {
          std::vector<uint32_t> new_deg(a.n);
          std::transform(coef_a.first.begin(), coef_a.first.end(),
                         coef_b.first.begin(), new_deg.begin(),
                         std::plus<uint32_t>());

          new_monomials.emplace(std::make_pair(new_deg, new_coef));
        }
      }
    }

    return PolynomialFF(a.n, new_monomials);
  }

  PolynomialFF PolynomialFF::mul_shift(const ff_map& a, const ff_map& b, uint32_t curr_deg) const {
    ff_map new_monomials;
    new_monomials.reserve(a.size()*b.size() + 1);

    for (const auto & coef_a : a) {
      for (const auto & coef_b : b) {

        FFInt new_coef = coef_a.second * coef_b.second;

        if (new_coef != 0) {
          std::vector<uint32_t> new_deg  = coef_a.first;
          new_deg[curr_deg] = coef_b.first[curr_deg];

          new_monomials.emplace(std::make_pair(new_deg, new_coef));
        }
      }
    }

    return PolynomialFF(n, new_monomials);
  }

  PolynomialFF PolynomialFF::add_shift(const std::vector<FFInt>& shift) const {
    if (shift.size() != n)
      throw std::runtime_error("Mismatch in sizes of the shift and variables!");

    PolynomialFF res;
    res.n = n;

    for (auto & mon : coefs) {
      PolynomialFF pow_poly;
      pow_poly.n = n;
      std::vector<uint32_t> powers = mon.first;

      for (uint32_t j = 0; j < n; ++j) {
        uint32_t deg = powers[j];
        FFInt tmp_shift = shift[j];

        // Calculate all terms originating from (x - a)^deg
        // by determining the binomial coefficients and adding
        // proper powers
        if (deg > 0) {
          ff_map tmp_pow_poly;

          if (tmp_shift > 0) {
            std::vector<std::vector<uint32_t>> tmp_powers(deg + 1, std::vector<uint32_t> (n));

            for (uint32_t k = 0; k <= deg; ++k) {
              tmp_powers[k][j] = deg - k;

              if (k == 0)
                tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], 1));
              else if (k == deg)
                tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], tmp_shift.pow(deg)));
              else
                tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], bin_coef(deg, k)*tmp_shift.pow(k)));
            }
          } else {
            std::vector<uint32_t> tmp_powers(n);
            tmp_powers[j] = deg;
            tmp_pow_poly.emplace(std::make_pair(tmp_powers, 1));
          }

          if (pow_poly.coefs.empty()) {
            pow_poly = PolynomialFF(n, tmp_pow_poly);
            pow_poly *= mon.second;
          } else
            pow_poly = mul_shift(pow_poly.coefs, tmp_pow_poly, j);
        }
      }

      if (!pow_poly.coefs.empty()) {
        //std::clock_t begin3 = clock();
        if (res.coefs.empty())
          res.coefs = pow_poly.coefs;
        else
          res += pow_poly;
      }
    }

    return res;
  }

  FFInt PolynomialFF::bin_coef(uint32_t n, uint32_t k) const {
    FFInt res = 1;

    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
      k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (uint32_t i = 0; i < k; ++i) {
      res *= (FFInt(n) - FFInt(i));
      res /= (FFInt(i) + FFInt(1));
    }

    return res;
  }

  /*void PolynomialFF::generate_hornerff() {
    if (!coefs.empty()) {
      std::vector<std::string> vars {};

      for (int i = 0; i < n; i++) {
        vars.emplace_back(std::string(1, ShuntingYardParser::get_var(i)));
      }

      s_y_fun.parse_function(generate_horner_coefs(0, coefs), vars);

      ff_map coefs_n_m_1 {};

      for (const auto & el : coefs) {
        std::vector<uint32_t> degs = el.first;
        degs.erase(degs.begin());

        if (coefs_n_m_1.find(degs) == coefs_n_m_1.end())
          coefs_n_m_1.emplace(std::make_pair(degs, el.second));
        else
          coefs_n_m_1[degs] += el.second;
      }

      s_y_fun.parse_function(generate_horner_coefs(0, coefs_n_m_1), vars);
    } else { // Interprete an empty polynomial as zero
      s_y_fun.parse_function("0", {"a"});
      s_y_fun_n_m_1.parse_function("0", {"a"});
      s_y_fun_map_n_m_1.parse_function("0", {"a"});
    }
  }

  std::string PolynomialFF::generate_horner_coefs(int index, const ff_map& monomials) {
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
        horner_coefs.emplace(std::make_pair(el.first, generate_horner_coefs(index + 1, el.second)));
      }

      uint32_t max_deg = tmp_coefs.begin() -> first;
      std::string horner_coef = "";

      if (max_deg > 0) {
        for (uint32_t i = 0; i < max_deg - 1; ++i) {
          horner_coef += "(";
        }

        const std::string var = std::string(1, ShuntingYardParser::get_var(index));

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

        const std::string var = std::string(1, ShuntingYardParser::get_var(index));

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
  }*/
}



