#include "Polynomial.hpp"

namespace firefly {

  Polynomial::Polynomial(std::vector<RationalNumber> coefs_) {
    for (uint i = 0; i < (uint) coefs_.size(); i++) {
      RationalNumber rn = coefs_.at(i);

      if (rn.numerator != 0) coefs.emplace_back(Monomial(rn, i));
    }
  }

  std::ostream &operator<<(std::ostream &out, const Polynomial &pol) {
    bool first = true;

    for (const auto mono : pol.coefs) {
      if (first) {
        out <<  mono.coef_rn << "*x^" << mono.deg;
        first = false;
      } else {
        out << " + " << mono.coef_rn << "*x^" << mono.deg;
      }

    }

    out << "\n";
    return out;
  }

}
