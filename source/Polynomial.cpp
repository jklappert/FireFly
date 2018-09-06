#include "Polynomial.hpp"

namespace firefly {

  Polynomial::Polynomial(const rn_map& coef) {
    for (const auto & el : coef) {
      coefs.emplace_back(Monomial(el.first, el.second));
    }
  }

  Polynomial::Polynomial(const Monomial& coef) {
    coefs.emplace_back(coef);
  }

  Polynomial::Polynomial() {}

  Polynomial Polynomial::operator*(const RationalNumber& rn) {
    for (auto & mon : coefs) {
      mon.coef = mon.coef * rn;
    }

    return *this;
  }

  Polynomial Polynomial::homogenize(uint degree) {
    for (auto & mon : coefs) {
      uint old_degree = 0;

      for (auto & power : mon.powers) old_degree += power;

      mon.powers.insert(mon.powers.begin(), degree - old_degree);
    }

    return *this;
  }

  Polynomial& Polynomial::operator-=(const Polynomial& b) {
    for (auto & coef_b : b.coefs) {
      int pos = -1;

      for (uint i = 0; i < (uint) coefs.size(); i++) {
        if (coef_b.powers == coefs[i].powers) pos = i;
      }

      if (pos == -1 && coef_b.coef.numerator != 0) {
        coefs.emplace_back(-coef_b);
      } else if (pos != -1) {
        coefs[pos].coef -= coef_b.coef;

        if (coefs[pos].coef.numerator == 0) coefs.erase(coefs.begin() + pos);
      }
    }

    return *this;
  }

  //todo can be optimized using an unordered_map
  Polynomial& Polynomial::operator+=(const Polynomial& b) {
    for (auto & coef_b : b.coefs) {
      int pos = -1;

      for (uint i = 0; i < (uint) coefs.size(); i++) {
        if (coef_b.powers == coefs[i].powers) pos = i;
      }

      if (pos == -1 && coef_b.coef.numerator != 0) {
        coefs.emplace_back(coef_b);
      } else if (pos != -1) {
        coefs[pos].coef += coef_b.coef;

        if (coefs[pos].coef.numerator == 0) coefs.erase(coefs.begin() + pos);
      }
    }

    return *this;
  }

  Polynomial& Polynomial::operator+=(const Monomial& b) {
    int pos = -1;

    for (uint i = 0; i < (uint) coefs.size(); i++) {
      if (b.powers == coefs[i].powers) pos = i;
    }

    if (pos == -1 && b.coef.numerator != 0) {
      coefs.emplace_back(b);
    } else if (pos != -1) {
      coefs[pos].coef += b.coef;

      if (coefs[pos].coef.numerator == 0) coefs.erase(coefs.begin() + pos);
    }

    return *this;
  }

  Polynomial Polynomial::operator*(const Polynomial& b) {
    rn_map new_monomials;
    new_monomials.reserve(coefs.size());

    for (auto & coef_a : coefs) {
      for (auto & coef_b : b.coefs) {
        Monomial new_coef = coef_a * coef_b;

        if (new_coef.coef.numerator != 0) {
          if (new_monomials.find(new_coef.powers) == new_monomials.end()) {
            new_monomials.emplace(std::make_pair(new_coef.powers, new_coef.coef));
          } else {
            new_monomials[new_coef.powers] += new_coef.coef;
          }
        }
      }
    }

    return Polynomial(new_monomials);
  }

  Polynomial Polynomial::operator*(const Monomial& b) {
    Polynomial a = *this;

    for (auto & coef_a : a.coefs) {
      coef_a = coef_a * b;
    }

    return a;
  }

  void Polynomial::sort() {
    std::sort(coefs.begin(), coefs.end());
  }

  void Polynomial::clear() {
    coefs.clear();
  }

  std::string Polynomial::string(const std::vector<std::string>& symbols) const {
    std::string str;

    for (const auto & mono : coefs) {
      str += mono.coef.string() + "*";

      for (uint i = 0; i < mono.powers.size(); i++) {
        if (mono.powers[i] > 1) {
          str += symbols[i] + "^" + std::to_string(mono.powers[i]) + "*";
        } else if (mono.powers[i] == 1) {
          str += symbols[i] + "*";
        }
      }

      str.erase(--str.end());
      str += "+";
    }

    str.erase(--str.end());
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

    out << "\n";
    return out;
  }

  PolynomialFF Polynomial::convert_to_PolynomialFF() {
    ff_map coefs_ff;
    uint n = coefs[0].powers.size();

    for (auto & coef : coefs) {
      mpz_class numerator = coef.coef.numerator % FFInt::p;

      if (numerator < 0) numerator += FFInt::p;

      mpz_class denominator = coef.coef.denominator % FFInt::p;
      FFInt coef_ff = FFInt(std::stoull(numerator.get_str())) / FFInt(std::stoull(denominator.get_str()));
      if(coef_ff.n > 0)
        coefs_ff.emplace(std::make_pair(coef.powers, coef_ff));
    }

    return PolynomialFF(n, coefs_ff);
  }

  Polynomial Polynomial::add_shift(std::vector<FFInt>& shift) {
    if (shift.size() != coefs[0].powers.size())
      throw std::runtime_error("Mismatch in sizes of the shift and variables!");

    uint n = shift.size();
    std::vector<RationalNumber> rn_shift(n);
    std::vector<uint> zero_deg(n);

    for (int i = 0; i < n; i++) {
      rn_shift[i] = RationalNumber(mpz_class(shift[i].n), 1);
    }

    Polynomial res;

    for (auto & mon : coefs) {
      Polynomial pow_poly;
      std::vector<uint> powers = mon.powers;
      std::vector<uint> decr_power = powers;

      for (int j = 0; j < n; j++) {
        uint deg = powers[j];

        if (deg > 0) {
          std::vector<uint> i_power(n);
          i_power[j] = 1;
          decr_power[j] = 0;
          rn_map add_shift;
          add_shift.emplace(std::make_pair(i_power, RationalNumber(1, 1)));
          add_shift.emplace(std::make_pair(zero_deg, rn_shift[j]));
          Polynomial tmp_pow_poly(add_shift);
          Polynomial mult_tmp_pow_poly = tmp_pow_poly;

          // TODO calc binomial coefficients to save some time
          for (int k = 1; k < deg; k++) {
            tmp_pow_poly = tmp_pow_poly * mult_tmp_pow_poly;
          }

          if (pow_poly.coefs.empty()) pow_poly = tmp_pow_poly;
          else pow_poly = pow_poly * tmp_pow_poly;
        }
      }
      // since always all variables are shifted decr_power := zero_deg
      if (!pow_poly.coefs.empty()) res += pow_poly * Monomial(decr_power, mon.coef);
    }

    return res;
  }
}

