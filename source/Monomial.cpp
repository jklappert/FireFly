#include "Monomial.hpp"

namespace firefly{
  Monomial::Monomial(RationalNumber &rn, uint deg_) : coef_rn(rn), deg(deg_) {
    is_coef_rn = true;
  }

  Monomial::Monomial(FFInt &ff, uint deg_) : coef_ff(ff), deg(deg_) {}
}