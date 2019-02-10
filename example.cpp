#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Tests.hpp"
#include "Logger.hpp"

using namespace firefly;

int main() {
  // Example for the reconstruction of a rational function
  uint n = 4;
  uint64_t prime = primes()[0];
  FFInt::set_new_prime(prime);
  RatReconst rec(n);
  //rec.set_tag("test");

  try {
    std::cout << "--------------------------------------------------------------\n";
    std::cout << "Interpolating rational function\n";
    std::cout << "--------------------------------------------------------------\n";
    // Initialize all values. t_yis are the scaled y-values and yis are the
    // unshifted ones
    std::vector<FFInt> t_yis(n);
    t_yis.reserve(n - 1);
    // t is the scaling variable
    FFInt t;

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    // Feed loop
    std::vector<FFInt> shift = rec.get_zi_shift_vec();

    while (!rec.is_done()) {
      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec.get_prime()) {
        if (!rec.need_shift()) {
          rec.disable_shift();
          shift = rec.get_zi_shift_vec();
          std::fill(shift.begin(), shift.end(), 0);
        }

        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec.get_prime());
        prime = primes()[rec.get_prime()];
        FFInt::set_new_prime(prime);
        rec.generate_anchor_points();
        kk = 0;
      }

      // Always set the scaling variable to a random value
      t = rec.get_rand();

      // Add the shift to the scaling variable
      FFInt z1 = t + shift[0];

      for (uint j = 2; j <= n; j++) {
        t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
      }

      std::vector<FFInt> yis(n);
      yis[0] = z1;

      for (uint j = 1; j < n; j++) {
        yis[j] = t_yis[j - 1];
      }

      //FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
      //FFInt num = n_eq_1(z1); // example for n = 1
      //FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
      FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients

      kk++;
      count++;

      // Feed the algorithm with the current zi_order
      rec.feed(t, num, rec.get_zi_order(), primes_used);
      rec.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec.get_result();
    std::cout << "--------------------------------------------------------------\n";
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  // Example for the reconstruction of a polynomial
  prime = primes()[0];
  FFInt::set_new_prime(prime);
  n = 3;
  PolyReconst rec_poly(n, 5);

  try {
    std::cout << "Interpolating polynomial\n";
    std::cout << "--------------------------------------------------------------\n";
    std::vector<FFInt> yis(n);
    rec_poly.generate_anchor_points();

    // Initialize some counters
    int count = 0;
    int kk = 0;
    uint primes_used = 0;

    // Reconstruction loop
    while (!rec_poly.is_done()) {

      // If a new prime is needed, set it, generate new random variables
      // and reset counters
      if (primes_used != rec_poly.get_prime()) {
        std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
        primes_used = std::max(primes_used, rec_poly.get_prime());

        prime = primes()[rec_poly.get_prime()];
        FFInt::set_new_prime(prime);
        rec_poly.generate_anchor_points();
        kk = 0;
      }

      yis = rec_poly.get_rand_zi_vec(rec_poly.get_zi_order());

      mpz_class cr_1_mpz;
      cr_1_mpz = "123456789109898799879870980";
      FFInt cr_1(cr_1_mpz);

      FFInt num = yis[0].pow(5) + yis[0] * yis[1].pow(4) + yis[0] * yis[1] * yis[2].pow(3) + yis[1].pow(5);

      kk++;
      count++;

      rec_poly.feed(num, rec_poly.get_zi_order(), primes_used);
      rec_poly.interpolate();
    }

    std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
    std::cout << rec_poly.get_result();
    std::cout << "--------------------------------------------------------------\n";
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  return 0;
}

