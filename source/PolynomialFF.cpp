#include "PolynomialFF.hpp"
#include <chrono>

namespace firefly {

  PolynomialFF::PolynomialFF(uint n_, ff_map coefs_) : n(n_), coefs(coefs_) {}

  PolynomialFF::PolynomialFF() {}

  FFInt PolynomialFF::calc(std::vector<FFInt> x) {
    FFInt res(0);

    for (const auto & term : coefs) {
      FFInt product(1);

      for (uint i = 0; i < n; i++) {
        product *= x[i].pow(term.first[i]);
      }

      res += term.second * product;
    }

    return res;
  }

  PolynomialFF PolynomialFF::operator+(const PolynomialFF& b) {
    PolynomialFF a = *this;
    ff_map new_coefs;

    if (a.coefs.size() > b.coefs.size()) {
      new_coefs = a.coefs;

      for (auto & el : b.coefs) {
        auto got = new_coefs.find(el.first);

        if (got == new_coefs.end()) {
          new_coefs.insert(el);
        } else {
          FFInt res = got -> second + el.second;

          if (res.n != 0) {
            got -> second = res;
          } else {
            new_coefs.erase(got -> first);
          }
        }
      }
    } else {
      new_coefs = b.coefs;

      for (auto & el : a.coefs) {
        auto got = new_coefs.find(el.first);

        if (got == new_coefs.end()) {
          new_coefs.insert(el);
        } else {
          FFInt res = got -> second + el.second;

          if (res.n != 0) {
            got -> second = res;
          } else {
            new_coefs.erase(got -> first);
          }
        }
      }
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::operator-(const PolynomialFF& b) {
    PolynomialFF a = *this;
    ff_map new_coefs;

    new_coefs = a.coefs;

    for (auto & el : b.coefs) {
      auto got = new_coefs.find(el.first);

      if (got == new_coefs.end()) {
        new_coefs.insert(std::make_pair(el.first, FFInt(0) - el.second));
      } else {
        FFInt res = got -> second - el.second;

        if (res.n != 0) {
          got -> second = res;
        } else {
          new_coefs.erase(got -> first);
        }
      }
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::operator*(const FFInt& ffint) {
    ff_map new_coefs;

    for (auto el : coefs) {
      new_coefs.insert(std::make_pair(el.first, el.second * ffint));
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::operator/(const FFInt& ffint) {
    ff_map new_coefs;

    for (auto el : coefs) {
      new_coefs.insert(std::make_pair(el.first, el.second / ffint));
    }

    return PolynomialFF(n, new_coefs);
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

  PolynomialFF PolynomialFF::mul(const uint zi) {
    ff_map new_coefs;

    for (auto coef_ : coefs) {
      std::vector<uint> new_element = coef_.first;
      new_element.at(zi - 1) ++;
      new_coefs.insert(std::make_pair(new_element, coef_.second));
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::homogenize(uint degree) {
    ff_map tmp_coef = coefs;
    coefs.clear();

    for (auto & mon : tmp_coef) {
      uint old_degree = 0;

      std::vector<uint> old_powers = mon.first;
      std::vector<uint> new_powers = old_powers;

      for (auto & power : old_powers) old_degree += power;

      new_powers.insert(new_powers.begin(), degree - old_degree);
      coefs[new_powers] = tmp_coef[old_powers];
    }

    tmp_coef.clear();
    return *this;
  }

  bool PolynomialFF::zero() {
    if (coefs.empty()) return true;

    //if(coef.size() == 1 && coef.begin()->second.n == 0) return true;
    return false;
  }

  std::vector<uint> PolynomialFF::max_deg() {
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

  std::vector<uint> PolynomialFF::min_deg() {
    max_deg();
    return min_degree;
  }

  PolynomialFF& PolynomialFF::operator-=(const PolynomialFF& b) {
    for (auto & coef_b : b.coefs) {
      if (coefs.find(coef_b.first) == coefs.end())
        coefs[coef_b.first] = FFInt(0) - coef_b.second;
      else
        coefs.at(coef_b.first) -= coef_b.second;
    }

    return *this;
  }

  //todo can be optimized using an unordered_map
  PolynomialFF& PolynomialFF::operator+=(const PolynomialFF& b) {
    for (auto & coef_b : b.coefs) {
      if (coefs.find(coef_b.first) == coefs.end())
        coefs[coef_b.first] = coef_b.second;
      else
        coefs.at(coef_b.first) += coef_b.second;
    }

    return *this;
  }

  PolynomialFF PolynomialFF::operator*(const PolynomialFF& b) {
    ff_map new_monomials;
    new_monomials.reserve(coefs.size());

    for (auto & coef_a : coefs) {
      for (auto & coef_b : b.coefs) {

        FFInt new_coef = coef_a.second * coef_b.second;

        if (new_coef != 0) {
          std::vector<uint> new_deg(n);
          std::transform(coef_a.first.begin(), coef_a.first.end(),
                         coef_b.first.begin(), new_deg.begin(),
                         std::plus<uint>());

          if (new_monomials.find(new_deg) == new_monomials.end()) {
            new_monomials.emplace(std::make_pair(new_deg, new_coef));
          } else {
            new_monomials[new_deg] += new_coef;
          }
        }
      }
    }

    return PolynomialFF(n, new_monomials);
  }

  PolynomialFF PolynomialFF::add_shift(std::vector<FFInt>& shift) {
    if (shift.size() != coefs.begin()->first.size())
      throw std::runtime_error("Mismatch in sizes of the shift and variables!");

    uint n = shift.size();
    std::vector<uint> zero_deg(n);

    PolynomialFF res;

    for (auto & mon : coefs) {
      PolynomialFF pow_poly;
      std::vector<uint> powers = mon.first;
      std::vector<uint> decr_power = powers;

//      std::clock_t begin = clock();

      for (uint j = 0; j < n; j++) {
        uint deg = powers[j];

        // Calculate all terms originating from (x - a)^deg
        // by determining the binomial coefficients and adding
        // proper powers
        if (deg > 0) {
          ff_map tmp_pow_poly;
          decr_power[j] = 0;
          std::vector<std::vector<uint>> tmp_powers(deg + 1, std::vector<uint> (n));

          for (uint jj = 0; jj <= deg; jj++) {
            tmp_powers[jj][j] = deg - jj;

            if (jj == 0) {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[jj], 1));
            } else if (jj == deg) {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[jj], shift[j].pow(deg)));
            } else {
              tmp_pow_poly.emplace(std::make_pair(tmp_powers[jj], bin_coef(deg, jj)*shift[j].pow(jj)));
            }
          }

          if (pow_poly.coefs.empty()) pow_poly = PolynomialFF(n, tmp_pow_poly);
          else pow_poly = pow_poly * PolynomialFF(n, tmp_pow_poly);
        }
      }

      //std::cout << " calculating terms took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";

      //begin = clock();

      // since always all variables are shifted decr_power := zero_deg
      if (!pow_poly.coefs.empty()) {
        ff_map monomial = {{decr_power, mon.second}};
        //pow_poly = pow_poly * PolynomialFF(n, monomial);
        //std::cout << " polynomial * monomial took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";
        //begin = clock();
        res += pow_poly * PolynomialFF(n, monomial);//* Monomial(decr_power, mon.second);
        //std::cout << " Poly + Poly took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";
      }
    }

    //std::cout << "size of poly " << res.coefs.size() << "\n";

    return res;
  }

  FFInt PolynomialFF::bin_coef(uint n, uint k) {
    FFInt res = 1;

    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
      k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (uint i = 0; i < k; ++i) {
      res *= (FFInt(n) - FFInt(i));
      res /= (FFInt(i) + FFInt(1));
    }

    return res;
  }

}



