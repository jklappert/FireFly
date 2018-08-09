#include "FFInt.hpp"
#include "iostream"

namespace firefly {

FFInt::FFInt(uint64_t n_, uint64_t p_) : n(n_), p(p_) {}

FFInt::FFInt(const FFInt& ffint) : n(ffint.n), p(ffint.p) {}

FFInt::FFInt() {}


FFInt& FFInt::operator+=(const FFInt& ffint){
   n += ffint.n;
   if(n >= p) n -= p;
   return *this;
}

FFInt& FFInt::operator-=(const FFInt& ffint){
   if(ffint.n > n) n += p;
   n -= ffint.n;
   return *this;
}

FFInt& FFInt::operator*=(const FFInt& ffint){
   n = mod_mul(n, ffint.n, p);
   return *this;
}

FFInt& FFInt::operator/=(const FFInt& ffint){
   n = mod_mul(n, mod_inv(ffint.n, p), p);
   return *this;
}

FFInt FFInt::pow(const FFInt& ffint) const{
   FFInt result;
   std::uint64_t exp;
   std::uint64_t base;

   if(ffint.n == 0){
      result.n = 1;
   } else{
      if(2*ffint.n < p){
	 // treat as positive exponent
	 exp = ffint.n;
	 base = n;
	 result.n = base;
      } else{
	 // treat as negative exponent
	 exp = p - ffint.n;
	 base = mod_inv(n, p); // =1/ffint1.n
	 result.n = base;
      }
      for(std::uint64_t i = 1; i != exp; ++i){
	result.n = mod_mul(result.n, base, p);
      }
   }
   return result;
}

FFInt FFInt::operator+(const FFInt& ffint){
   auto sum = ffint.n + n;
   if(sum >= p) sum -= p;
   return FFInt(sum, p);
}

FFInt FFInt::operator-(const FFInt& ffint){
   auto diff = n;
   if(ffint.n > diff) diff += p;
   diff -= ffint.n;
   return FFInt(diff, p);
}

FFInt FFInt::operator*(const FFInt& ffint){
   return FFInt(mod_mul(n, ffint.n, p), p);
}

FFInt FFInt::operator/(const FFInt& ffint){
   return FFInt(mod_mul(n, mod_inv(ffint.n, p), p), p);
}

bool FFInt::operator==(const FFInt& ffint){
   return (n == ffint.n);
}

bool FFInt::operator!=(const FFInt& ffint){
   return (n != ffint.n);
}

FFInt& FFInt::operator=(const FFInt& ffint){
   n = ffint.n;
   p = ffint.p;
   return *this;
}

uint64_t FFInt::mod_mul(uint64_t a, uint64_t b, const uint64_t p) const{
   long double x;
   uint64_t c;
   int64_t r;
   if(a >= p) a %= p;
   if(b >= p) b %= p;
   x = a;
   c = x*b/p;
   r = (int64_t)(a*b - c*p) % (int64_t)p;
   return r < 0 ? r + p : r;
}

uint64_t FFInt::mod_inv(const uint64_t a, const uint64_t p) const{
   int64_t t{0};
   int64_t newt{1};
   int64_t tmpt;
   uint64_t r{p};
   uint64_t newr{a};
   uint64_t tmpr;
   uint64_t q;
   while (newr){
      q = r/newr;
      tmpt = t;
      t = newt;
      newt = tmpt - q*newt;
      tmpr = r;
      r = newr;
      newr = tmpr - q*newr;
   }
   // if(r > 1) throw init_error("mod_inv: not invertible");
   return t < 0 ? t + p : t;
}

}