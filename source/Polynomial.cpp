#include "Polynomial.hpp"

namespace firefly {
   
Polynomial::Polynomial(std::vector<FFInt> coef_) : coef(coef_){
   deg = coef.size() - 1;
}

Polynomial::Polynomial() {}


Polynomial Polynomial::operator+(const Polynomial& b){
   Polynomial a = *this;
   std::vector<FFInt> newCoefs;
   if(a.deg >= b.deg){
      newCoefs = a.coef;
      for(int i = 0; i <= b.deg; i++) newCoefs.at(i) += b.coef.at(i);
   } else{
      newCoefs = b.coef;
      for(int i = 0; i <= a.deg; i++) newCoefs.at(i) += a.coef.at(i);
   }
   return Polynomial(newCoefs);
}

Polynomial Polynomial::operator-(const Polynomial& b){
   Polynomial a = *this;
   std::vector<FFInt> newCoefs;
   if(a.deg >= b.deg){
      newCoefs = a.coef;
      for(int i = 0; i <= b.deg; i++) newCoefs.at(i) -= b.coef.at(i);
   } else{
      newCoefs = b.coef;
      FFInt zero (0, a.coef.at(1).p);
      for(int i = 0; i <= a.deg; i++) newCoefs.at(i) = a.coef.at(i) - newCoefs.at(i);
      for(int i = a.deg + 1; i <= b.deg; i++) newCoefs.at(i) = zero - newCoefs.at(i);
   }
   return Polynomial(newCoefs);
}

Polynomial Polynomial::operator*(const Polynomial& b){
   Polynomial a = *this;
   std::vector<FFInt> newCoefs;
   const double newDeg = a.deg + b.deg;
   for(int i = 0; i <= newDeg; i++){
      FFInt ffint(0, a.coef.at(0).p);
      for(int j = 0; j <= i; j++){
	 if(a.deg >= j && b.deg >= (i - j)){ 
	    ffint += a.coef.at(j)*b.coef.at(i - j);
	 }
      }
      newCoefs.push_back(ffint);
   }
   return Polynomial(newCoefs);
}

Polynomial& Polynomial::operator=(const Polynomial& a){
   coef = a.coef;
   deg = a.deg;
   return *this;
}

std::ostream& operator<<(std::ostream& out, const Polynomial& a){
   for(int i = 0; i < (int) a.coef.size(); i++){
      const uint64_t n = a.coef.at(i).n;
      if(n != 0) out << " + " << n << "*x^" << i;
   }
   return out;
}

}