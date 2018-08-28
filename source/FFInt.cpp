#include <sstream>
#include "FFInt.hpp"

namespace firefly {

  uint64_t FFInt::p;

  FFInt::FFInt(const uint64_t n_) : n(n_) {}

  FFInt::FFInt(const FFInt& ffint) : n(ffint.n) {}

  FFInt::FFInt(const std::string& str, const std::vector<std::pair<std::string, uint64_t>>& replacements) {
    for (const auto & var : replacements) {
      if (var.first == str) {
        n = var.second;
        return;
      }
    }

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

  FFInt::FFInt() {}

  FFInt& FFInt::operator+=(const FFInt& ffint) {
    n += ffint.n;

    if (n >= p) n -= p;

    return *this;
  }

  FFInt& FFInt::operator-=(const FFInt& ffint) {
    if (ffint.n > n) n += p;

    n -= ffint.n;
    return *this;
  }

  FFInt& FFInt::operator*=(const FFInt& ffint) {
    n = mod_mul(n, ffint.n);
    return *this;
  }

  FFInt& FFInt::operator/=(const FFInt& ffint) {
    n = mod_mul(n, mod_inv(ffint.n));
    return *this;
  }

  FFInt FFInt::pow(const FFInt& ffint) const {
    FFInt result;
    std::uint64_t exp;
    std::uint64_t base;

    if (ffint.n == 0) {
      result.n = 1;
    } else {
      if (2 * ffint.n < p) {
        // treat as positive exponent
        exp = ffint.n;
        base = n;
        result.n = base;
      } else {
        // treat as negative exponent
        exp = p - ffint.n;
        base = mod_inv(n);  // =1/ffint1.n
        result.n = base;
      }

      for (std::uint64_t i = 1; i != exp; ++i) {
        result.n = mod_mul(result.n, base);
      }
    }

    return result;
  }

  FFInt FFInt::operator+(const FFInt& ffint) {
    auto sum = ffint.n + n;

    if (sum >= p) sum -= p;

    return FFInt(sum);
  }

  FFInt FFInt::operator-(const FFInt& ffint) {
    auto diff = n;

    if (ffint.n > diff) diff += p;

    diff -= ffint.n;
    return FFInt(diff);
  }

  FFInt FFInt::operator-() {
    return FFInt(p - n);
  }

  /*FFInt FFInt::operator*(const FFInt& ffint) {
    return FFInt(mod_mul(n, ffint.n));
  }*/

  FFInt FFInt::operator/(const FFInt& ffint) {
    return FFInt(mod_mul(n, mod_inv(ffint.n)));
  }

  bool FFInt::operator==(const FFInt& ffint) const {
    return (n == ffint.n);
  }

  bool FFInt::operator!=(const FFInt& ffint) const {
    return (n != ffint.n);
  }

  uint64_t FFInt::mod_mul(uint64_t a, uint64_t b) const {
    long double x;
    uint64_t c;
    int64_t r;

    if (a >= p) a %= p;

    if (b >= p) b %= p;

    x = a;
    c = x * b / p;
    r = (int64_t)(a * b - c * p) % (int64_t) p;
    return r < 0 ? r + p : r;
  }

  uint64_t FFInt::mod_inv(const uint64_t a) const {
    int64_t t {0};
    int64_t newt {1};
    int64_t tmpt;
    uint64_t r {p};
    uint64_t newr {a};
    uint64_t tmpr;
    uint64_t q;

    while (newr) {
      q = r / newr;
      tmpt = t;
      t = newt;
      newt = tmpt - q * newt;
      tmpr = r;
      r = newr;
      newr = tmpr - q * newr;
    }

    // if(r > 1) throw init_error("mod_inv: not invertible");
    return t < 0 ? t + p : t;
  }

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
      if (result) result = mod_mul(result, 1000000000000000000uLL);

      result += intchunk;
      result %= p;
    }

    return result;
  }

  FFInt operator*(const FFInt& a, const FFInt& b) {
    return FFInt(a.mod_mul(a.n, b.n));
  }

  FFInt pow(const FFInt& ffint, const FFInt& power) {
    return ffint.pow(power);
  }

  std::ostream& operator<<(std::ostream& out, const FFInt& ffint) {
    out << ffint.n;
    return out;
  }

}
