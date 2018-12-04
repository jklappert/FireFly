#include "Polynomial.hpp"
#include <chrono>

namespace firefly {

  Polynomial::Polynomial(const rn_map& coef) {
    coefs = coef;
    n = coefs.begin()->first.size();
  }

  Polynomial::Polynomial(const Monomial& coef) {
    coefs[coef.powers] = coef.coef;
    n = coef.powers.size();
  }

  Polynomial::Polynomial() {}

  Polynomial Polynomial::operator*(const RationalNumber& rn) {
    for (auto & mon : coefs) {
      mon.second = mon.second * rn;
    }

    return *this;
  }

  Polynomial Polynomial::homogenize(uint degree) {
    rn_map tmp_coef = coefs;
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

  Polynomial& Polynomial::operator-=(const Polynomial& b) {
    for (auto & coef_b : b.coefs) {
      try {
        coefs.at(coef_b.first) -= coef_b.second;
      } catch (std::out_of_range& e) {
        coefs[coef_b.first] = RationalNumber(0, 1);
        coefs[coef_b.first] -= coef_b.second;
      }
    }

    return *this;
  }

  //todo can be optimized using an unordered_map
  Polynomial& Polynomial::operator+=(const Polynomial& b) {
    for (auto & coef_b : b.coefs) {
      try {
        coefs.at(coef_b.first) += coef_b.second;
      } catch (std::out_of_range& e) {
        coefs[coef_b.first] = RationalNumber(0, 1);
        coefs[coef_b.first] += coef_b.second;
      }
    }

    return *this;
  }

  Polynomial& Polynomial::operator+=(const Monomial& b) {
    try {
      coefs.at(b.powers) += b.coef;
    } catch (std::out_of_range& e) {
      coefs[b.powers] = RationalNumber(0, 1);
      coefs[b.powers] += b.coef;
    }

    return *this;
  }

  Polynomial Polynomial::operator*(const Polynomial& b) {
    rn_map new_monomials;
    new_monomials.reserve(coefs.size());

    for (auto & coef_a : coefs) {
      for (auto & coef_b : b.coefs) {
        Monomial new_coef = Monomial(coef_a.first, coef_a.second) * Monomial(coef_b.first, coef_b.second);

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
    Polynomial c;

    for (auto & coef_a : a.coefs) {
      Monomial a_mon(coef_a.first, coef_a.second);
      a_mon = a_mon * b;
      c.coefs[a_mon.powers] = a_mon.coef;
    }

    a.clear();
    return c;
  }

  void Polynomial::sort() {
    //std::sort(coefs.begin(), coefs.end());
  }

  void Polynomial::clear() {
    coefs.clear();
  }

  std::string Polynomial::string(const std::vector<std::string>& symbols) const {
    std::string str;

    for (const auto & mono : coefs) {
      str += mono.second.string() + "*";

      for (uint i = 0; i < mono.first.size(); i++) {
        if (mono.first[i] > 1) {
          str += symbols[i] + "^" + std::to_string(mono.first[i]) + "*";
        } else if (mono.first[i] == 1) {
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
        out <<  mono.second << "*x^(";

        for (const auto i : mono.first) {
          out << i << ",";
        }

        out << "\b)";
        first = false;
      } else {
        out << " + " << mono.second << "*x^(";

        for (const auto i : mono.first) {
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
    uint n = coefs.begin()->first.size();

    for (auto & coef : coefs) {
      mpz_class numerator = coef.second.numerator % FFInt::p;

      if (numerator < 0) numerator += FFInt::p;

      mpz_class denominator = coef.second.denominator % FFInt::p;
      FFInt coef_ff = FFInt(std::stoull(numerator.get_str())) / FFInt(std::stoull(denominator.get_str()));

      if (coef_ff.n > 0)
        coefs_ff.emplace(std::make_pair(coef.first, coef_ff));
    }

    return PolynomialFF(n, coefs_ff);
  }

  Polynomial Polynomial::add_shift(std::vector<FFInt>& shift) {
    if (shift.size() != coefs.begin()->first.size())
      throw std::runtime_error("Mismatch in sizes of the shift and variables!");

    uint n = shift.size();
    std::vector<RationalNumber> rn_shift(n);
    std::vector<uint> zero_deg(n);

    for (uint i = 0; i < n; i++) {
      rn_shift[i] = RationalNumber(mpz_class(shift[i].n), 1);
    }

    Polynomial res;

    for (auto & mon : coefs) {
      Polynomial pow_poly;
      std::vector<uint> powers = mon.first;
      std::vector<uint> decr_power = powers;

      std::clock_t begin = clock();

      for (uint j = 0; j < n; j++) {
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
          for (uint k = 1; k < deg; k++) {
            tmp_pow_poly = tmp_pow_poly * mult_tmp_pow_poly;
          }

          if (pow_poly.coefs.empty()) pow_poly = tmp_pow_poly;
          else pow_poly = pow_poly * tmp_pow_poly;
        }
      }
      
      std::cout << "pow Poly!\n" << pow_poly;

      std::cout << " calculating terms took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";

      begin = clock();

      // since always all variables are shifted decr_power := zero_deg

      if (!pow_poly.coefs.empty()) {
        pow_poly = pow_poly * Monomial(decr_power, mon.second);
        std::cout << " polynomial * monomial took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";
        begin = clock();
        res += pow_poly;//* Monomial(decr_power, mon.second);
        std::cout << " Poly + Poly took : " << float(clock() - begin) / CLOCKS_PER_SEC << "\n";
      }
    }

    std::cout << "size of poly " << res.coefs.size() << "\n";

    return res;
  }
}

