#include <cstdlib>
#include <iostream>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"

namespace firefly {
   
PolyReconst::PolyReconst(int n_) : n(n_) {}

std::vector<FFInt> PolyReconst::reconst(){
   const int64_t prime = primes().at(0);
   int maxDegree = 1000;
   const int breakCondition = 3;

   yi.emplace_back(FFInt(std::rand() % prime, prime));
   ai.emplace_back(num(prime, yi.back()));
   
   for(int i = 1; i < maxDegree; i++){
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      ai.emplace_back(compAi(i, i, num(prime, yi.back())));
      if(ai.at(i).n == 0) {
	 if(i > breakCondition){
	    bool nonZero = false;
	    for(int j = ai.size(); j > ai.size() - breakCondition; j--){
	       if(ai.at(j - 1).n != 0){
		  nonZero = true;
		  break;
	       }
	    }
	    if(!nonZero) break;
	 }
      }
      if(i == maxDegree - 1) maxDegree += 1000;
   }
   return ai;
}

FFInt PolyReconst::compAi(int i, int ip, const FFInt& num){
   if(ip == 0){
      return num;
   } else{
      return (compAi(i, ip - 1, num) - ai.at(ip - 1))/(yi.at(i) - yi.at(ip - 1));
   }
}

FFInt PolyReconst::num(uint64_t p, const FFInt& y){
   FFInt a0 (3, p);
   FFInt a1 (6, p);
   FFInt a2 (18, p);
   FFInt exp (2, p);
   return a0 + a1*y + a2*y.pow(exp);
}

}