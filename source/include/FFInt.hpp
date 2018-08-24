#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace firefly {

  class FFInt {
  public:
    /**
     *    A constructor
     *    @param n_ an integer which is a member of the finite fild
     *    @param p_ a prime number
     */
    FFInt(const uint64_t n_);
    /**
     *    A constructor
     *    @param ffint a FFInt object
     */
    FFInt(const FFInt& ffint);
    FFInt(const std::string& str, const std::vector<std::pair<std::string, uint64_t>>& replacements);
    /**
     *    Default constructor
     */
    FFInt();

    // defining new operators for finite field arithmetic
    FFInt& operator=(const FFInt&) = default;
    FFInt& operator+=(const FFInt&);
    FFInt& operator-=(const FFInt&);
    FFInt& operator*=(const FFInt&);
    FFInt& operator/=(const FFInt&);
    FFInt operator+(const FFInt&);
    FFInt operator-(const FFInt&);
    FFInt operator-();
    FFInt operator*(const FFInt&);
    FFInt operator/(const FFInt&);
    bool operator==(const FFInt&) const;
    bool operator!=(const FFInt&) const;
    FFInt pow(const FFInt& ffint) const;

    uint64_t n; /**< the integer member of the finite field */
    static uint64_t p; /**< the prime defining the finite field */
  private:
    /**
     *    A function to calculate a*b mod p
     *    taken from https://en.wikipedia.org/wiki/Modular_arithmetic
     *    @param a the input
     *    @param b the ceil
     *    @param p the prime which defines the finite field
     */
    uint64_t mod_mul(uint64_t a, uint64_t b) const;
    /**
     *    Extended Euclidian algorithm to calculate the multiplicative
     *    invrse.
     *    mod_inv(a,p) solves a*t = 1 mod p for t.
     */
    uint64_t mod_inv(const uint64_t a) const;
    uint64_t parse_longint(const std::string& str);
  };

  FFInt operator*(const FFInt& a, const FFInt& b);
  FFInt pow(const FFInt& ffint, const FFInt& power);
  std::ostream& operator<<(std::ostream& out, const FFInt& ffint);

}
