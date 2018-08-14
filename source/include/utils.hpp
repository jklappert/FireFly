#pragma once
#include <gmpxx.h>
#include "RationalNumber.hpp"

namespace firefly{
  std::pair<mpz_class, mpz_class> chineseRemainder (
    const std::pair<mpz_class, mpz_class> &p1,
    const std::pair<mpz_class, mpz_class> &p2);

  RationalNumber getRationalCoef (const mpz_class& a, const mpz_class& p);
}