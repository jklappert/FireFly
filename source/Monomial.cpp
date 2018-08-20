#include "Monomial.hpp"

namespace firefly{

  Monomial::Monomial(const std::vector<uint> &powers_, const RationalNumber &coef_) : powers(powers_), coef(coef_) {}

  bool Monomial::operator<(const Monomial& m2) {
    uint deg1 = 0;
    uint deg2 = 0;
    
    for (uint i = 0; i < powers.size(); i++) {
      deg1 += powers[i];
      deg2 += m2.powers[i];
    }
    
    if (deg1 < deg2) {
      return true;
    } else {
      return deg1 == deg2 && powers < m2.powers;
    }
  }

  bool Monomial::operator>(const Monomial& m2) {
    uint deg1 = 0;
    uint deg2 = 0;
    
    for (uint i = 0; i < powers.size(); i++) {
      deg1 += powers[i];
      deg2 += m2.powers[i];
    }
    
    if (deg1 > deg2) {
      return true;
    } else {
      return deg1 == deg2 && powers > m2.powers;
    }
  }
  
}