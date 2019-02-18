#include "PolynomialFF.hpp"
#include <chrono>

namespace firefly {

  PolynomialFF::PolynomialFF(uint32_t n_, ff_map coefs_) : n(n_), coefs(coefs_) {}

  PolynomialFF::PolynomialFF() {}

  FFInt PolynomialFF::calc(std::vector<FFInt> x) const {
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

  FFInt PolynomialFF::calc_n_m_1(std::vector<FFInt> x) const {
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
      if (coefs.find(coef_b.first) == coefs.end())
        coefs[coef_b.first] = FFInt(0) - coef_b.second;
      else
        coefs[coef_b.first] -= coef_b.second;
    }

    return *this;
  }

  PolynomialFF& PolynomialFF::operator+=(const PolynomialFF& b) {
    for (const auto & coef_b : b.coefs) {
      if (coefs.find(coef_b.first) == coefs.end())
        coefs[coef_b.first] = coef_b.second;
      else
        coefs[coef_b.first] += coef_b.second;
    }

    return *this;
  }

  PolynomialFF operator*(const PolynomialFF& a, const FFInt& ffint) {
    ff_map new_coefs = a.coefs;

    for (auto& el : new_coefs) {
      el.second *= ffint;
    }

    return PolynomialFF(a.n, new_coefs);
  }


  PolynomialFF operator/(const PolynomialFF& a, const FFInt& ffint) {
    ff_map new_coefs = a.coefs;
    FFInt inv = 1 / ffint;

    for (auto& el : new_coefs) {
      el.second *= inv;
    }

    return PolynomialFF(a.n, new_coefs);
  }

  PolynomialFF& PolynomialFF::operator*=(const FFInt& ffint) {
    for (auto& el : coefs) {
      el.second *= ffint;
    }

    return *this;
  }

  PolynomialFF& PolynomialFF::operator/=(const FFInt& ffint) {
    FFInt inv = 1 / ffint;

    for (auto& el : coefs) {
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

  PolynomialFF PolynomialFF::mul(const uint32_t zi) {
    ff_map new_coefs {};

    for (const auto& coef_ : coefs) {
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

  bool PolynomialFF::zero() {
    if (coefs.empty())
      return true;
    else if(coefs.size() == 1 && coefs.begin()->second == 0)
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
    new_monomials.reserve(a.coefs.size() + b.coefs.size());

    for (const auto & coef_a : a.coefs) {
      for (const auto & coef_b : b.coefs) {

        FFInt new_coef = coef_a.second * coef_b.second;

        if (new_coef != 0) {
          std::vector<uint32_t> new_deg(a.n);
          std::transform(coef_a.first.begin(), coef_a.first.end(),
                         coef_b.first.begin(), new_deg.begin(),
                         std::plus<uint32_t>());

          if (new_monomials.find(new_deg) == new_monomials.end())
            new_monomials.emplace(std::make_pair(new_deg, new_coef));
          else
            new_monomials[new_deg] += new_coef;
        }
      }
    }

    return PolynomialFF(a.n, new_monomials);
  }


  PolynomialFF PolynomialFF::add_shift(const std::vector<FFInt>& shift) {
    if (shift.size() != n)
      throw std::runtime_error("Mismatch in sizes of the shift and variables!");

    PolynomialFF res;
    res.n = n;
    //std::clock_t begin2 = clock();

    for (auto & mon : coefs) {
      PolynomialFF pow_poly;
      std::vector<uint32_t> powers = mon.first;

      for (uint32_t j = 0; j < n; ++j) {
        uint32_t deg = powers[j];

        // Calculate all terms originating from (x - a)^deg
        // by determining the binomial coefficients and adding
        // proper powers
        if (deg > 0) {
          ff_map tmp_pow_poly;
          std::vector<std::vector<uint32_t>> tmp_powers(deg + 1, std::vector<uint32_t> (n));

          for (uint32_t k = 0; k <= deg; ++k) {
            tmp_powers[k][j] = deg - k;

            if (k == 0) {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], 1));
            } else if (k == deg) {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], shift[j].pow(deg)));
            } else {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[k], bin_coef(deg, k)*shift[j].pow(k)));
            }
          }

          if (pow_poly.coefs.empty())
            pow_poly = PolynomialFF(n, tmp_pow_poly);
          else
            pow_poly = pow_poly * PolynomialFF(n, tmp_pow_poly);
        }
      }

      // since always all variables are shifted decr_power := zero_deg
      if (!pow_poly.coefs.empty()){
        pow_poly *= mon.second;
        res += pow_poly;
      }
    }
        //std::cout << " shift took : " << float(clock() - begin2) / CLOCKS_PER_SEC << "\n";

    return res;
  }

  FFInt PolynomialFF::bin_coef(uint32_t n, uint32_t k) {
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

}


