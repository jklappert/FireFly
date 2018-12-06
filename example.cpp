#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
//#include "utils.hpp"
#include <algorithm>
//#include <ctime>
//#include <Eigen/Sparse>

int main() {
  const uint n = 4;
  //  typedef Eigen::Matrix<firefly::FFInt, n, 1> vec;
  //  std::vector<std::vector<firefly::FFInt>> coef_mat;
  uint64_t prime = firefly::primes()[0];
  firefly::FFInt::p = prime;
  firefly::RatReconst rec(n);

  // construct data
  /*  Eigen::SparseMatrix<firefly::FFInt> mat (n,n);
  vec b, x;

  std::clock_t begin = clock();
  std::cout << "Writing equations...\n";
  for(int i = 0; i < n; i++){
    std::vector<firefly::FFInt> eq (n + 1, 0);
    for (int j = i; j < n; j++){
      firefly::FFInt rand1 = firefly::FFInt(std::rand());
      //      mat.insert(i,j) = rand1.n;
      eq[j] = rand1;
    }
    firefly::FFInt rand2 = firefly::FFInt(std::rand());
    //    b[i] = rand2.n;
    eq[n] = rand2;
    coef_mat.emplace_back(eq);
  }

  std::cout << "Finished writing equations in " << float(clock() - begin) / CLOCKS_PER_SEC << " s.\n";

  begin = clock();
  std::cout << "Solving Gauss system...\n";
  std::vector<firefly::FFInt> res = firefly::solve_gauss_system(coef_mat.size(), coef_mat);
  std::cout << "Solved Gauss system in " << float(clock() - begin) / CLOCKS_PER_SEC << " s.\nResults:\n";
  /*for(const auto & el : res){
    std::cout << el << "\n";
  }*/
  /*std::cout << "End of results.\n";
  Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
  solver.analyzePattern(mat);
  //  Eigen::ConjugateGradient<Eigen::SparseMatrix<firefly::FFInt>, Eigen::Lower|Eigen::Upper> cg;
  begin = clock();
  std::cout << "Solving sparse system with Eigen...\n";
  //  x = cg.compute(mat).solve(b);
  std::cout << "Solved sparse system with Eigen in " << float(clock() - begin) / CLOCKS_PER_SEC << " s.\nResults:\n";
  /*for(int i = 0; i < n; i++){
    std::cout << x[i] << "\n";
  }*/
  //  std::cout << "End of results.\n";

  // Example for the reconstruction of a rational function
  try {
    // Initialize all values. t_yis are the scaled y-values and yis are the
    // unshifted ones
    std::vector<firefly::FFInt> t_yis;
    std::vector<firefly::FFInt> yis;
    yis.reserve(n - 1);
    t_yis.reserve(n - 1);
    // t is the scaling variable
    firefly::FFInt t;
    // Get and set a prime number
    prime = firefly::primes()[rec.get_prime()];
    firefly::FFInt::p = prime;
    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    // Feed loop
    while (!rec.is_done()) {
      if (yis.size() > 0) {
        for (uint j = 2; j <= n; j++) {
          t_yis[j - 2] = yis[j - 2];
        }
      }

      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec.get_prime()) {
        rec.disable_shift();
        rec.generate_anchor_points();

        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec.get_prime());

        prime = firefly::primes()[rec.get_prime()];
        firefly::FFInt::p = prime;

        kk = 0;

        for (uint j = 2; j <= n; j++) {
          yis[j - 2] = rec.rand_zi[std::make_pair(j, rec.get_zi_order()[j - 2])];
          t_yis[j - 2] = yis[j - 2];
        }
      }

      // Always set the scaling variable to a random value
      t = rec.get_rand();

      // Fill the t_yis vector if it is empty
      if (t_yis.size() == 0) {
        for (uint j = 2; j <= n; j++) {
          yis.emplace_back(rec.rand_zi[std::make_pair(j, rec.get_zi_order()[j - 2])]);

          if (primes_used == 0)
            t_yis.emplace_back(t * yis[j - 2] + firefly::RatReconst::shift[j - 1]);
          else
            t_yis.emplace_back(yis[j - 2]);
        }
      }

      // Check whether we are currently reconstructing multivariate polynomials,
      // change the yis and t_yis accordingly and add the static shift
      if (rec.get_zi() >= 2) {
        for (uint j = 2; j <= n; j++) {
          yis[j - 2] = rec.rand_zi[std::make_pair(j, rec.get_zi_order()[j - 2])];

          if (primes_used == 0)
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          else
            t_yis[j - 2] = yis[j - 2];
        }
      } else {
        for (uint j = 2; j <= n; j++) {
          if (primes_used == 0)
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          else
            t_yis[j - 2] = yis[j - 2];
        }
      }

      // Add the shift to the scaling variable
      firefly::FFInt z1 = t + firefly::RatReconst::shift[0];

      // Some examples for number for which one needs to use the Chinese
      // Remainder theorem
      mpz_class cr_1_mpz;
      cr_1_mpz = "123456789109898799879870980";
      mpz_class cr_2_mpz;
      cr_2_mpz = "123456789109898799879";
      firefly::FFInt cr_1(cr_1_mpz);
      firefly::FFInt cr_2(cr_2_mpz);

      // example for n = 4 using the Chinese Remainder theorem
      firefly::FFInt den = cr_1 * (((z1.pow(3) - 12 * z1.pow(2) + 48 * z1 - 64) * t_yis[1].pow(2))
                                   * t_yis[0].pow(5) + ((-3 * z1.pow(3) + 36 * z1.pow(2)
                                                         - 144 * z1 + 192) * t_yis[1].pow(2)) * t_yis[0].pow(4) + ((2 * z1.pow(3) - 24 * z1.pow(2)
                                                             + 96 * z1 - 128) * t_yis[1].pow(2)) * t_yis[0].pow(3) + ((2 * z1.pow(3) - 24 * z1.pow(2)
                                                                 + 96 * z1 - 128) * t_yis[1].pow(2)) * t_yis[0].pow(2) + ((-3 * z1.pow(3) + 36 * z1.pow(2)
                                                                     - 144 * z1 + 192) * t_yis[1].pow(2)) * t_yis[0] + (z1.pow(3) - 12 * z1.pow(2) + 48 * z1 - 64)
                                   * t_yis[1].pow(2));
//       firefly::FFInt den = 1;
      firefly::FFInt num = cr_2 * ((-6 * z1.pow(3) + 54 * z1.pow(2) - 156 * z1 + 144)
                                   * t_yis[0].pow(4) + ((-4 * z1.pow(3) + 36 * z1.pow(2) - 104 * z1 + 96) * t_yis[1]
                                                        + 9 * z1.pow(3) - 84 * z1.pow(2) + 252 * z1 - 240) * t_yis[0].pow(3) + ((46 * z1.pow(3)
                                                            - 389 * z1.pow(2) + 1074 * z1 - 960) * t_yis[1] - 3 * z1.pow(3) + 30 * z1.pow(2) - 96 * z1 + 96)
                                   * t_yis[0].pow(2) + ((-10 * z1.pow(3) + 93 * z1.pow(2) - 278 * z1 + 264)
                                                        * t_yis[1]) * t_yis[0]) + 0*z1.pow(10)*t_yis[0].pow(10)*t_yis[1].pow(10)*t_yis[2].pow(10);
      /*firefly::FFInt num = t_yis[0].pow(100) + t_yis[1].pow(100);
      for(int i = 1; i < 100; i++){
        num += t_yis[0].pow(i)*t_yis[1].pow(100-i);
      }*/
      /*std::clock_t begin = clock();
      for(int k = 1; k < 50; k++){
        for(kk = k; kk > 0; kk --){
          num += t_yis[0].pow(k)*t_yis[1].pow(kk);
        }
      }*/

      //firefly::FFInt num = z1.pow(4) + 3*t_yis[0].pow(5) + t_yis[1].pow(2);
      //firefly::FFInt den = 2*z1*t_yis[0]*t_yis[1].pow(2) + 3*t_yis[0];

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
        tmp_vec = rec.get_zi_order();

      // Feed the algorithm with the current zi_order
      rec.feed(t, num / den, tmp_vec, primes_used);
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec.get_result();
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  return 0;
}

