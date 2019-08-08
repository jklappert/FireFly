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
//=================================================================================
#include "FFThieleInterpolator.hpp"

namespace firefly {

  ThieleInterpolator::ThieleInterpolator() {}

  bool ThieleInterpolator::add_point(const FFInt& num, const FFInt& yi) {
    bool check = false;
    ti.emplace_back(yi);
    const uint32_t i = ti.size() - 1;

    if (i == 0) {
      ai.emplace_back(num);
    } else {
      if (num == comp_fyi(i - 1, ti.back())) check = true;

      if (!check)
        ai.emplace_back(comp_ai(i, num));
    }

    return check;
  }

  std::pair<ff_map, ff_map> ThieleInterpolator::get_result() {
    ti.pop_back();
    return construct_canonical();
  }

  FFInt ThieleInterpolator::comp_ai(uint32_t i, const FFInt& num) {
    FFInt res = num;

    if (i > 0) {
      for (uint32_t ip_tmp = 1; ip_tmp != i + 1; ip_tmp++) {
        res = (ti[i] - ti[ip_tmp - 1]) / (res - ai[ip_tmp - 1]);
      }
    }

    return res;
  }

  FFInt ThieleInterpolator::comp_fyi(uint32_t i, const FFInt& y) {
    FFInt res = ai[i];

    if (i > 0) {
      for (uint32_t ip_tmp = 1; ip_tmp != i + 1; ip_tmp++) {
        res = ai[i - ip_tmp] + (-ti[i - ip_tmp] + y) / res;
      }
    }

    return res;
  }

  std::pair<ff_map, ff_map> ThieleInterpolator::construct_canonical() {
    if (ai.size() == 1) {
      ff_map numerator_ff;
      std::vector<uint32_t> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[0]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(numerator_ff, denominator_ff);
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical();
      FFInt mti = -ti[0];
      return std::make_pair((r.first * ai[0] + r.second * mti + r.second.mul(1)).coefs,
                            r.first.coefs);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> ThieleInterpolator::iterate_canonical() {
    ff_map numerator_ff;
    std::vector<uint32_t> zero_deg = {0};
    numerator_ff.emplace(std::make_pair(zero_deg, ai.back()));
    ff_map denominator_ff;
    denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
    std::pair<PolynomialFF, PolynomialFF> res = std::make_pair(PolynomialFF(1, numerator_ff),
                                                               PolynomialFF(1, denominator_ff));

    for (uint32_t i_tmp = ai.size() - 2; i_tmp != 0; i_tmp--) {
      FFInt mti = -ti[i_tmp];
      res = std::pair<PolynomialFF, PolynomialFF> (res.first * ai[i_tmp] + res.second.mul(1) + res.second * mti,
                                                   res.first);
    }

    return res;
  }

  ThieleInterpolator& ThieleInterpolator::operator=(const ThieleInterpolator& other) {
    if (this != &other) {
      ai = other.ai;
      ti = other.ti;
    }

    return *this;
  }

  ThieleInterpolator& ThieleInterpolator::operator=(ThieleInterpolator && other) {
    if (this != &other) {
      ai = std::move(other.ai);
      ti = std::move(other.ti);
    }

    return *this;
  }

  ThieleInterpolator::ThieleInterpolator(const ThieleInterpolator& other) {
    ai = other.ai;
    ti = other.ti;
  }

  ThieleInterpolator::ThieleInterpolator(ThieleInterpolator && other) {
    ai = std::move(other.ai);
    ti = std::move(other.ti);
  }
}
