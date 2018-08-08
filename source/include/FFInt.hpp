#pragma once

#include <cstdint>

namespace firefly{
   
   class FFInt{
   public:
      /**
       * 	A constructor
       * 	@param n_ an integer which is a member of the finite fild
       * 	@param p_ a prime number
       */
      FFInt(std::uint64_t n_, std::uint64_t p_);
      /**
       * 	A constructor
       * 	@param ffint a FFInt object
       */
      FFInt(const FFInt& ffint);
      FFInt();
      
      // defining new operators for finite field arithmetic
      FFInt& operator =(const FFInt&);
      FFInt& operator +=(const FFInt&);
      FFInt& operator -=(const FFInt&);
      FFInt& operator *=(const FFInt&);
      FFInt& operator /=(const FFInt&);
      FFInt operator +(const FFInt&);
      FFInt operator -(const FFInt&);
      FFInt operator *(const FFInt&);
      FFInt operator /(const FFInt&);
      FFInt pow(const FFInt&);
      bool operator ==(const FFInt&);
      bool operator !=(const FFInt&);
      
      uint64_t n; /**< the integer member of the finite field */
      uint64_t p; /**< the prime defining the finite field */
   private:
      /**
       * 	A function to calculate a*b mod p
       * 	taken from https://en.wikipedia.org/wiki/Modular_arithmetic
       * 	@param a the input
       * 	@param b the ceil
       * 	@param p the prime which defines the finite field
       */
      uint64_t mod_mul(uint64_t a, uint64_t b, const uint64_t p);
      /**
       * 	Extended Euclidian algorithm to calculate the multiplicative
       * 	invrse.
       * 	mod_inv(a,p) solves a*t = 1 mod p for t.
       */
      uint64_t mod_inv(const uint64_t a, const uint64_t p);
   };
}