#include "PolynomialFF.hpp"

namespace firefly {

  PolynomialFF::PolynomialFF(uint n_, ff_map coef_) : n(n_), coef(coef_) {}

  PolynomialFF::PolynomialFF() {}

  FFInt PolynomialFF::calc(std::vector<FFInt> x) {
    FFInt res(0);

    for (auto & term : coef) {
      FFInt product(1);

      for (uint i = 0; i < x.size(); i++) {
        product *= x[i].pow(FFInt(term.first[i]));
      }

      res += term.second * product;
    }

    return res;
  }

  PolynomialFF PolynomialFF::operator+(const PolynomialFF& b) {
    PolynomialFF a = *this;
    ff_map new_coefs;

    if (a.coef.size() > b.coef.size()) {
      new_coefs = a.coef;

      for (auto & el : b.coef) {
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
      new_coefs = b.coef;

      for (auto & el : a.coef) {
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

    new_coefs = a.coef;

    for (auto & el : b.coef) {
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

    for (auto el : coef) {
      new_coefs.insert(std::make_pair(el.first, el.second * ffint));
    }

    return PolynomialFF(n, new_coefs);
  }

  PolynomialFF PolynomialFF::operator/(const FFInt& ffint) {
    ff_map new_coefs;

    for (auto el : coef) {
      new_coefs.insert(std::make_pair(el.first, el.second / ffint));
    }

    return PolynomialFF(n, new_coefs);
  }

  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a) {
    for (auto & coef_ : a.coef) {
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

    for (auto coef_ : coef) {
      std::vector<uint> new_element = coef_.first;
      new_element.at(zi - 1) ++;
      new_coefs.insert(std::make_pair(new_element, coef_.second));
    }

    return PolynomialFF(n, new_coefs);
  }


  bool PolynomialFF::zero() {
    if (coef.empty()) return true;

    //if(coef.size() == 1 && coef.begin()->second.n == 0) return true;
    return false;
  }

  std::vector<uint> PolynomialFF::max_deg() {
    if (max_degree.empty()) {
      int tmp_max;
      int tmp_min;
      for (const auto c : coef) {
        int tmp_deg = 0;

        for (const auto i : c.first) {
          tmp_deg = + i;
        }

        if (max_degree.empty()) {
          tmp_max = tmp_deg;
          tmp_min = tmp_deg;
        }

        tmp_max = std::max(tmp_deg, tmp_max);
        tmp_min = std::min(tmp_deg, tmp_min);

        if(tmp_max == tmp_deg) max_degree = c.first;
        if(tmp_min == tmp_deg) min_degree = c.first;

      }
    }
    return max_degree;
  }

  std::vector<uint> PolynomialFF::min_deg() {
    max_deg();
    return min_degree;
  }
}
