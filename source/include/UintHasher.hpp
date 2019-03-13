//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

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
