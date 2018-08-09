#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"

namespace firefly {

RatReconst::RatReconst(int n_) : n(n_), prime(primes().at(0)) {}

std::vector<FFInt> RatReconst::reconst(){
   uint maxDegree = 1000;
   const int breakCondition = 3;

   yi.emplace_back(FFInt(std::rand() % prime, prime));
   ai.emplace_back(num(prime, yi.back()));
   fyi.emplace_back(ai.at(0));

   for(uint i = 1; i < maxDegree; i++){
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      fyi.emplace_back(num(prime, yi.back()));
      ai.emplace_back(compAi(i, i, fyi.back()));
      if(fyi.back() == compFyi(i - 1, yi.back())){
	 if(i > breakCondition){
	    bool nonequal = false;
	    for(uint j = ai.size(); j > ai.size() - breakCondition; j--){
	       if(fyi.at(j - 1) != compFyi(j - 2, yi.at(j - 1))){
		  nonequal = true;
		  break;
	       }
	    }
	    if(!nonequal) break;
	 }
      }
      if(i == maxDegree - 1) maxDegree += 1000;
   }
   for(int i = 0; i < breakCondition; i++){
      yi.pop_back();
      ai.pop_back();
      fyi.pop_back();
   }
   return ai;
}

FFInt RatReconst::compAi(int i, int ip, const FFInt& num){
   if(ip == 0){
      return num;
   } else{
      return (yi.at(i) - yi.at(ip - 1))/(compAi(i, ip - 1, num) - ai.at(ip - 1));
   }
}

FFInt RatReconst::compFyi(int i, FFInt& y){
   return iterateCanonicalNum(i, i, y);
}

FFInt RatReconst::iterateCanonicalNum(uint i, uint ip, const FFInt& y){
   if(ip == 0){
      return ai.at(i);
   } else{
      return ai.at(i - ip) + (FFInt(0, prime) - yi.at(i - ip) + y)/iterateCanonicalNum(i, ip - 1, y);
   }
}

void RatReconst::constrCanonical(){
   /*if(ai.size() == 0){
      INFO_MSG("Polynomial not yet reconstructed or 0.");
   } else if(ai.size() == 1){
      std::vector<FFInt> coef{ai.at(0)};
      Polynomial poly(coef);
      canonical = poly;
   } else{
      std::vector<FFInt> coef{ai.at(0)};
      Polynomial poly(coef);
      canonical = poly + iterateCanonical(1);
   }*/
}

Polynomial RatReconst::iterateCanonical(uint i){
   /*if(i < ai.size() - 1){
      std::vector<FFInt> coef1{ (FFInt(0, prime) - yi.at(i - 1))/ai.at(i),
	 FFInt(1, prime)/ai.at(i) };
      std::vector<FFInt> coef2{ FFInt(0, prime) - yi.at(i - 1), FFInt(1, prime) };
      Polynomial poly1(coef1);
      Polynomial poly2(coef2);
      return poly1 + poly2 / iterateCanonical(i + 1);
   } else{
      std::vector<FFInt> coef1{ (FFInt(0, prime) - yi.at(i - 1))/ai.at(i),
	 FFInt(1, prime)/ai.at(i) };
      Polynomial poly1(coef1);
      return poly1;
   }*/
   return Polynomial();
}


FFInt RatReconst::num(uint64_t p, const FFInt& y){
   FFInt a0 (3, p);
   FFInt a1 (6, p);
   FFInt a2 (18, p);
   FFInt a3 (25, p);
   FFInt a4 (30, p);
   FFInt a5 (2, p);
   FFInt a6 (7, p);
   FFInt a7 (100, p);
   FFInt a8 (13, p);
   FFInt exp2 (2, p);
   FFInt exp3 (3, p);
   FFInt exp4 (4, p);
   FFInt exp5 (5, p);
   FFInt exp6 (6, p);
   FFInt exp7 (7, p);
   FFInt exp8 (8, p);

   return (a0 + a1*y + a4*y.pow(exp2) + a4*y.pow(exp3))/(a2 + a3*y );
}

}