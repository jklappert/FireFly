#include "Polynomial.hpp"

namespace firefly {

  Polynomial::Polynomial(const rn_map &coef) {
    for (const auto & el : coef) {
      coefs.emplace_back(Monomial(el.first, el.second));
    }

    std::sort(coefs.begin(), coefs.end());
  }

  Polynomial::Polynomial() {}

  Polynomial Polynomial::operator*(const RationalNumber& rn) {
    for(auto& mon : coefs){
      mon.coef = mon.coef * rn;
    }
    return *this;
  }

  std::ostream &operator<<(std::ostream &out, const Polynomial &pol) {
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

}
