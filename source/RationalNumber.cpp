#include "RationalNumber.hpp"

namespace firefly {
  
  RationalNumber::RationalNumber (mpz_class nominator_, mpz_class denominator_) : nominator(nominator_), denominator(denominator_) {}
  
  std::ostream &operator<< (std::ostream &out, const RationalNumber &a) {
    out << "(" << a.nominator.get_str() << "/" << a.denominator.get_str() << ")";
    return out;
  }


}