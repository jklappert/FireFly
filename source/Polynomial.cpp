#include "Polynomial.hpp"

namespace firefly {

  Polynomial::Polynomial(const rn_map& coef) {
    for (const auto & el : coef) {
      coefs.emplace_back(Monomial(el.first, el.second));
    }
  }

  Polynomial::Polynomial() {}

  Polynomial Polynomial::operator*(const RationalNumber& rn) {
    for (auto & mon : coefs) {
      mon.coef = mon.coef * rn;
    }

    return *this;
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

  Polynomial Polynomial::homogenize(uint degree) {
    for (auto & mon : coefs) {
      uint old_degree = 0;

      for (auto & power : mon.powers) old_degree += power;

      mon.powers.insert(mon.powers.begin(), degree - old_degree);
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
      } else {
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
    } else {
      coefs[pos].coef += b.coef;

      if (coefs[pos].coef.numerator == 0) coefs.erase(coefs.begin() + pos);
    }

    return *this;
  }

  Polynomial Polynomial::operator*(const Polynomial& b) {
    //std::vector<Monomial> new_monomials;
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
}

