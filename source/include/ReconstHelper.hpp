// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <cstdint>
#include <vector>
#include <string>

namespace firefly {

  /**
   *    A collection of 63-bit primes
   */
  const std::vector<uint64_t>& primes();

  /**
   * A collection of 127-bit primes
   */
  const std::vector<std::string>& primes_128();
}
