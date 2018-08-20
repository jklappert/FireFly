#pragma once

#include <vector>

namespace firefly {
  class UintHasher {
  public:
    std::size_t operator()(std::vector<uint> const &vec) const {
      std::size_t seed = vec.size();

      for (auto & i : vec) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }

      return seed;
    }
  };
}
