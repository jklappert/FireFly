// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================
#include <sstream>
#include "FFInt.hpp"
#ifdef FLINT
#include "flint/ulong_extras.h"
#endif

namespace firefly {

  uint64_t FFInt::p;
  uint64_t FFInt::p_inv;

  FFInt::FFInt(const FFInt& ffint) : n(ffint.n) {}

  FFInt::FFInt(mpz_class in) {
    in = in % p;

    if (in < 0) in = p + in;

    n = std::stoull(in.get_str());
  }

  // copied from Kira
  FFInt::FFInt(const std::string& str, const std::vector<std::pair<std::string, uint64_t>>& replacements) {
    for (const auto & var : replacements) {
      if (var.first == str) {
        if (var.second > p) {
          n = var.second % p;
        } else {
          n = var.second;
        }

        return;
      }
    }

    if (str.front() == '-') throw std::runtime_error("Negative number\n");

    std::istringstream ss {str};
    auto success = static_cast<bool>(ss >> n);

    if (!(success && ss.rdbuf()->in_avail() == 0)) {
      if (str.empty()) {
        // (ss >> n) fails if ss is empty
        throw std::runtime_error("FFInt: empty argument");
      } else if (std::isalpha(str.front())) {
        throw std::runtime_error("Unkown or invalid coefficient string \"" + str + "\"");
      } else {
        n = parse_longint(str);
      }
    } else if (n >= p) {
      // special case: the parsed value fits into n, but is >= p
      n %= p;
    }
  }

  FFInt::FFInt() {
    n = 0;
  }

#ifdef FLINT
  FFInt& FFInt::operator+=(const FFInt& ffint) {
    n = n_addmod(n, ffint.n, p);

    return *this;
  }

  FFInt& FFInt::operator-=(const FFInt& ffint) {
    n = n_submod(n, ffint.n, p);

    return *this;
  }

  FFInt& FFInt::operator*=(const FFInt& ffint) {
    n = n_mulmod2_preinv(n, ffint.n, p, p_inv);
    return *this;
  }

  FFInt& FFInt::operator/=(const FFInt& ffint) {
    n = n_mulmod2_preinv(n, n_invmod(ffint.n, p), p, p_inv);
    return *this;
  }

  FFInt FFInt::pow(const FFInt& ffint) const {
    return FFInt(n_powmod2_preinv(n, ffint.n, p, p_inv));
  }

  FFInt operator+(const FFInt& a, const FFInt& b) {
    return FFInt(n_addmod(a.n, b.n, FFInt::p));
  }

  FFInt operator-(const FFInt& a, const FFInt& b) {
    return FFInt(n_submod(a.n, b.n, FFInt::p));
  }
#endif

  FFInt FFInt::operator-() const {
    return FFInt(p - n);
  }

  FFInt FFInt::operator+() const {
    return FFInt(n);
  }

  FFInt FFInt::operator++() {
    FFInt tmp(n);
    tmp += 1;
    return tmp;
  }

  FFInt FFInt::operator++(int) {
    FFInt tmp(n);
    tmp += 1;
    return tmp;
  }

  FFInt FFInt::operator--(int) {
    FFInt tmp(n);
    tmp -= 1;
    return tmp;
  }

  FFInt FFInt::operator--() {
    FFInt tmp(n);
    tmp -= 1;
    return tmp;
  }

  bool FFInt::operator!() const {
    return !n;
  }

  // copied from Kira
  uint64_t FFInt::parse_longint(const std::string& str) {
    // Parse a long integer, passed as a string, take the modulus wrt. prime
    // and return it. The string is split into chunks of at most 18 digits
    // which are then put together via modular arithmetic.
    // Return zero if the string is empty.
    //
    // Make sure the input is an unsigned integer without whitespace padding.
    for (const auto ch : str) {
      if (!std::isdigit(ch)) throw std::runtime_error("parse_longint(): invalid number string \"" + str + "\"");
    }

    uint64_t result = 0;
    std::size_t pos = 0;
    std::size_t len = ((str.size() - 1) % 18) + 1;

    while (pos < str.size()) {
      std::string strchunk = str.substr(pos, len);
      pos += len;
      len = 18;
      uint64_t intchunk;
      std::istringstream ss {strchunk};
      ss >> intchunk;

      // result=0 in the first pass or when the string is zero padded
      // on the left so that the first (few) chunks give zero.
      if (result) result =  (FFInt(result) * FFInt(1000000000000000000uLL)).n;//n_mulmod2_preinv(result, 1000000000000000000uLL, FFInt::p, FFInt::p_inv);

      result = (FFInt(result) + FFInt(intchunk)).n;// n_addmod(result, intchunk, FFInt::p);
    }

    return result;
  }

  bool operator==(const FFInt& a, const FFInt& b) {
    return (a.n == b.n);
  }

  bool operator!=(const FFInt& a, const FFInt& b) {
    return (a.n != b.n);
  }

  bool operator<(const FFInt& a, const FFInt& b) {
    return (a.n < b.n);
  }

  bool operator<=(const FFInt& a, const FFInt& b) {
    return (a.n <= b.n);
  }

  bool operator>(const FFInt& a, const FFInt& b) {
    return (a.n > b.n);
  }

  bool operator>=(const FFInt& a, const FFInt& b) {
    return (a.n >= b.n);
  }

#ifdef FLINT
  FFInt operator/(const FFInt& a, const FFInt& b) {
    return FFInt(n_mulmod2_preinv(a.n, n_invmod(b.n, FFInt::p), FFInt::p, FFInt::p_inv));
  }

  FFInt operator*(const FFInt& a, const FFInt& b) {
    return FFInt(n_mulmod2_preinv(a.n, b.n, FFInt::p, FFInt::p_inv));
  }
#endif

  FFInt pow(const FFInt& ffint, const FFInt& power) {
    return ffint.pow(power);
  }

  std::ostream& operator<<(std::ostream& out, const FFInt& ffint) {
    out << ffint.n;
    return out;
  }
#ifdef FLINT
  void FFInt::set_new_prime(uint64_t prime) {
    FFInt::p = prime;
    FFInt::p_inv = n_preinvert_limb(prime);
  }
#else
  void FFInt::set_new_prime(uint64_t prime) {
    FFInt::p = prime;
  }
#endif

  void firefly_exists(void) {}
}
