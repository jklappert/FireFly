#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <gmpxx.h>

namespace firefly {

  class FFInt {
  public:
    /**
     *    A constructor
     *    @param n_ an integer which is a member of the finite fild
     *    @param p_ a prime number
     */
    template<typename T>
    FFInt(const T n_);
    /**
     *    A constructor
     *    @param ffint a FFInt object
     */
    FFInt(const FFInt& ffint);
    FFInt(mpz_class& in);
    FFInt(const std::string& str, const std::vector<std::pair<std::string, uint64_t>>& replacements);
    /**
     *    Default constructor
     */
    FFInt();
    static void set_new_prime(uint64_t prime);

    // defining new operators for finite field arithmetic
    FFInt& operator=(const FFInt&) = default;
    FFInt& operator+=(const FFInt&);
    FFInt& operator-=(const FFInt&);
    FFInt& operator*=(const FFInt&);
    FFInt& operator/=(const FFInt&);
    FFInt operator+(const FFInt&);
    FFInt operator-(const FFInt&);
    FFInt operator-();
    //FFInt operator*(const FFInt&);
    FFInt operator/(const FFInt&);
    bool operator==(const FFInt&) const;
    bool operator!=(const FFInt&) const;
    FFInt pow(const FFInt& ffint) const;

    uint64_t n; /**< the integer member of the finite field */
    static uint64_t p; /**< the prime defining the finite field */
    static uint64_t p_inv;
  private:
    uint64_t parse_longint(const std::string& str);
  };

  FFInt operator*(const FFInt&, const FFInt&);
  FFInt pow(const FFInt& ffint, const FFInt& power);
  std::ostream& operator<<(std::ostream& out, const FFInt& ffint);

  template<typename T>
  FFInt::FFInt(const T n_) {
    if (n_ >= 0) {
      if ((uint64_t) n_ < p) {
        n = n_;
      } else {
        n = n_ % p;
      }
    } else if (n_ < 0) {
      n = p - (-n_) % p;
    }
  }
}
