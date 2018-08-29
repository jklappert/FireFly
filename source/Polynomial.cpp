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
    for(auto & mon : coefs){
      uint old_degree = 0;
      for(auto & power : mon.powers) old_degree += power;
      mon.powers.insert(mon.powers.begin(), degree - old_degree);
    }
    return *this;
  }

  //todo can be optimized using an unordered_map
  Polynomial& Polynomial::operator+=(const Polynomial& a) {
    for(auto & coef_a : a.coefs){
      int pos = -1;
      for(uint i = 0; i < (uint) coefs.size(); i++){
        if(coef_a.powers == coefs[i].powers) pos = i;
      }
      if(pos == -1) {
        coefs.emplace_back(coef_a);
      } else{
        coefs[pos].coef += coef_a.coef;
      }
    }
    return *this;
  }

  void Polynomial::sort() {
    std::sort(coefs.begin(), coefs.end());
  }

}
