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
      if (num == comp_fyi(i - 1, i - 1, ti.back())) check = true;

      if (!check)
        ai.emplace_back(comp_ai(i, i, num));
    }

    return check;
  }

  std::pair<ff_map, ff_map> ThieleInterpolator::get_result() {
    ti.pop_back();
    return construct_canonical();
  }

  FFInt ThieleInterpolator::comp_ai(int i, int ip, const FFInt& num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt ai_i = comp_ai(i, ip - 1, num);
      return (ti[i] - ti[ip - 1]) / (ai_i - ai[ip - 1]);
    }
  }

  FFInt ThieleInterpolator::comp_fyi(uint32_t i, uint32_t ip, const FFInt& y) {
    if (ip == 0) {
      return ai[i];
    } else {
      return ai[i - ip] + (-ti[i - ip] + y) / comp_fyi(i, ip - 1, y);
    }
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
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(1);
      FFInt mti = -ti[0];
      return std::make_pair((r.first * ai[0] + r.second * mti + r.second.mul(1)).coefs,
                            r.first.coefs);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> ThieleInterpolator::iterate_canonical(uint32_t i) {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(i + 1);
      FFInt mti = -ti[i];
      return std::pair<PolynomialFF, PolynomialFF> (fnp1.first * ai[i] + fnp1.second.mul(1) + fnp1.second * mti,
                                                    fnp1.first);
    } else {
      ff_map numerator_ff;
      std::vector<uint32_t> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[i]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    }
  }
}
