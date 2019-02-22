// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <vector>

namespace firefly {
  class UintHasher {
  public:
    std::size_t operator()(std::vector<uint32_t> const& vec) const {
      std::size_t seed = vec.size();

      for (auto & i : vec) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }

      return seed;
    }
  };

  class UintPairHasher {
  public:
    std::size_t operator()(std::pair<uint32_t, uint32_t> const& pair) const {
      std::size_t seed = 2;

      seed ^= pair.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= pair.second + 0x9e3779b9 + (seed << 6) + (seed >> 2);

      return seed;
    }
  };
}
