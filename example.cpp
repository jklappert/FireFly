#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include <algorithm>

int main() {
  uint n = 4;
  uint64_t prime = firefly::primes()[0];
  firefly::FFInt::p = prime;
  firefly::RatReconst rec_1(n);
  std::unordered_map<int, firefly::RatReconst> map;

  for (int m = 0; m < 1; m++) {
    firefly::RatReconst tmp_rec(n);
    map.emplace(m, std::move(tmp_rec));
  }

  try {

    for (int m = 0; m < 1; m++) {
      std::vector<firefly::FFInt> t_yis;
      std::vector<firefly::FFInt> yis;
      yis.reserve(n - 1);
      t_yis.reserve(n - 1);
      firefly::FFInt t;
      firefly::RatReconst& rec = map.at(m);
      prime = firefly::primes()[rec.prime_number];
      firefly::FFInt::p = prime;
      int count = 0;
      int kk = 0;
      uint primes_used = 0;

      while (!rec.done) {
        if (yis.size() > 0) {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = yis[j - 2];
          }
        }

        if (primes_used != rec.prime_number) {
          rec.generate_anchor_points();

          std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
          primes_used = std::max(primes_used, rec.prime_number);
          prime = firefly::primes()[rec.prime_number];
          firefly::FFInt::p = prime;

          kk = 0;

          for (uint j = 2; j <= n; j++) {
            yis[j - 2] = rec.rand_zi[std::make_pair(j, rec.curr_zi_order[j - 2])];
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        }

        t = rec.get_rand();

        if (t_yis.size() == 0) {
          for (uint j = 2; j <= n; j++) {
            yis.emplace_back(rec.rand_zi[std::make_pair(j, rec.curr_zi_order[j - 2])]);
            t_yis.emplace_back(t * yis[j - 2] + firefly::RatReconst::shift[j - 1]);
          }
        }

        if (rec.zi >= 2) {
          for (uint j = 2; j <= n; j++) {
            yis[j - 2] = rec.rand_zi[std::make_pair(j, rec.curr_zi_order[j - 2])];
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        } else {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        }

        firefly::FFInt z1 = t + firefly::RatReconst::shift[0];

        // Some examples for number for which one needs to use the Chinese
        // Remainder theorem
        mpz_class cr_1_mpz;
        cr_1_mpz = "123456789109898799879870980";
        mpz_class cr_2_mpz;
        cr_2_mpz = "123456789109898799879";
        firefly::FFInt cr_1(cr_1_mpz);
        firefly::FFInt cr_2(cr_2_mpz);

        // example for n = 4
        firefly::FFInt den = cr_1*(((z1.pow(3)-12*z1.pow(2)+48*z1-64)*t_yis[1].pow(2))
          *t_yis[0].pow(5)+((-3*z1.pow(3)+36*z1.pow(2)
          -144*z1+192)*t_yis[1].pow(2))*t_yis[0].pow(4)+((2*z1.pow(3)-24*z1.pow(2)
          +96*z1-128)*t_yis[1].pow(2))*t_yis[0].pow(3)+((2*z1.pow(3)-24*z1.pow(2)
          +96*z1-128)*t_yis[1].pow(2))*t_yis[0].pow(2)+((-3*z1.pow(3)+36*z1.pow(2)
          -144*z1+192)*t_yis[1].pow(2))*t_yis[0]+(z1.pow(3)-12*z1.pow(2)+48*z1-64)
          *t_yis[1].pow(2));
        firefly::FFInt num = cr_2*((-6*z1.pow(3)+54*z1.pow(2)-156*z1+144)
          *t_yis[0].pow(4)+((-4*z1.pow(3)+36*z1.pow(2)-104*z1+96)*t_yis[1]
          +9*z1.pow(3)-84*z1.pow(2)+252*z1-240)*t_yis[0].pow(3)+((46*z1.pow(3)
          -389*z1.pow(2)+1074*z1-960)*t_yis[1]-3*z1.pow(3)+30*z1.pow(2)-96*z1+96)
          *t_yis[0].pow(2)+((-10*z1.pow(3)+93*z1.pow(2)-278*z1+264)
          *t_yis[1])*t_yis[0]);

        // example for n = 1
        /*firefly::FFInt num = (576 * z1.pow(12) - 35145 * z1.pow(11)
                              + 946716 * z1.pow(10) - 14842335 * z1.pow(9)
                              + 150236238 * z1.pow(8) - 1028892363 * z1.pow(7)
                              + 4853217576 * z1.pow(6) - 15724949577 * z1.pow(5)
                              + 34208917206 * z1.pow(4) - 47506433412 * z1.pow(3)
                              + 37933483608 * z1.pow(2) - 13296184128 * z1 + 71850240);
        firefly::FFInt den = (16 * z1.pow(12) - 960 * z1.pow(11)
                              + 25456 * z1.pow(10) - 393440 * z1.pow(9)
                              + 3934768 * z1.pow(8) - 26714240 * z1.pow(7)
                              + 125545488 * z1.pow(6) - 408157280 * z1.pow(5)
                              + 899198016 * z1.pow(4) - 1278172800 * z1.pow(3)
                              + 1055033856 * z1.pow(2) - 383201280 * z1);*/
        kk++;
        count++;
        std::vector<uint> tmp_vec;

        if (n > 1)
          tmp_vec = std::vector<uint>(rec.curr_zi_order.begin(), rec.curr_zi_order.end() - 1);

        rec.feed(t, num / den, tmp_vec);
      }

      std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
      std::cout << rec.get_result();
    }
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  return 0;
}
