#include "Monomial.hpp"

namespace firefly {

  Monomial::Monomial(const std::vector<uint>& powers_, const RationalNumber& coef_) : powers(powers_), coef(coef_) {}

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
      // since (0,0,1,0) > (0,1,0,0) we need to change the relational operator
      return deg1 == deg2 && powers > m2.powers;
    }
  }

  bool Monomial::operator>(const Monomial& m2) {
    uint deg1 = 0;
    uint deg2 = 0;

    for (uint i = 0; i < (uint) powers.size(); i++) {
      deg1 += powers[i];
      deg2 += m2.powers[i];
    }

    if (deg1 > deg2) {
      return true;
    } else {
      // since (0,0,1,0) > (0,1,0,0) we need to change the relational operator
      return deg1 == deg2 && powers < m2.powers;
    }
  }

  Monomial Monomial::operator*(const Monomial& b) {
    Monomial a = *this;
    a.coef = a.coef*b.coef;
    for(uint i = 0; i < (uint) powers.size(); i++){
      a.powers[i] += b.powers[i];
    }
    return a;
  }

  Monomial Monomial::operator-() const {
    return Monomial(powers, - coef);
  }
}
