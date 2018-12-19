#pragma once

#include <cstdint>
#include <vector>
#include <string>

namespace firefly {

  /**
   *    A collection of 63 bit prime numbers
   */
  const std::vector<uint64_t>& primes();
  
  /**
   * 
   */
  const std::vector<std::string>& primes_128();
}
