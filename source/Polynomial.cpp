#include "Polynomial.hpp"

namespace firefly {

  Polynomial::Polynomial(std::vector<RationalNumber> coefs_) : coefs(coefs_) {}

  std::ostream &operator<<(std::ostream &out, const Polynomial &pol) {
    for (int i = 0; i < (int) pol.coefs.size() ; i++) {
      const RationalNumber ri = pol.coefs.at(i);
      const mpz_class ci = ri.numerator;

      if (i == 0 && ri.numerator != 0) {
        out << ri;
      } else if (ci != 0) {
        out << " + " << ri << "*x^" << i;
      }
    }

    out << "\n";
    return out;
  }

}
