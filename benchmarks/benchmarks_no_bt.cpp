//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "firefly/Reconstructor.hpp"

namespace firefly {
  class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
  public:
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
      std::vector<FFIntTemp> result = par.evaluate_pre(values);
      return result;
    }

    inline void prime_changed() {
      par.precompute_tokens();
    }

  private:
    ShuntingYardParser par;
  };
}

using namespace firefly;
// Performs benchmarks of hep-ph:2004.01463
int main() {
  INFO_MSG("Performing benchmarks of hep-ph:2004.01463");
  PolyReconst::use_bt = false;
  std::string root_dir = FIREFLY_ROOT_DIR;
  // Eq. (18-20) & Eq. (29) from hep-ph:2004.01463
  INFO_MSG("Using hybrid racer without Ben-Or/Tiwari");
  for (int i = 1; i != 5; ++i) {
    switch (i) {
      case 1: {
        INFO_MSG("Eq. (18)");
        std::cout << "-----------------------------------------------\n";
        break;
      }
      case 2: {
        INFO_MSG("Eq. (19)");
        std::cout << "-----------------------------------------------\n";
	break;
      }
      case 3: {
	INFO_MSG("Eq. (20)");
        std::cout << "-----------------------------------------------\n";
	break;
      }
      case 4: {
	INFO_MSG("Eq. (29)");
        std::cout << "-----------------------------------------------\n";
	break;
      }
    }
    if (i == 1) {
      ShuntingYardParser par(root_dir + "/benchmarks/f" + std::to_string(i) + ".m", {"x1", "x2", "x3", "x4",
            "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20"});
      BlackBoxUser bb(par);
      Reconstructor<BlackBoxUser> reconst(20, 1, 1, bb);
      reconst.enable_shift_scan();
      reconst.reconstruct();
    } else {
      ShuntingYardParser par(root_dir + "/benchmarks/f" + std::to_string(i) + ".m", {"x1", "x2", "x3", "x4", "x5"});
      BlackBoxUser bb(par);
      Reconstructor<BlackBoxUser> reconst(5, 1, 1, bb);
      reconst.enable_shift_scan();
      reconst.reconstruct();
    }
    RatReconst::reset();
    std::cout << "-----------------------------------------------\n";
  }
  return 0;
}
