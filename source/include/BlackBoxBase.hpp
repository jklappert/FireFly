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

namespace firefly {
  /**
   * @class BlackBoxBase
   * @brief The base class of the black box
   */
  template<typename BlackBoxTemp>
  class BlackBoxBase {
  public:
    /**
     * TODO templates
     *  The evaluation of the black box. This function is called from Reconstructor.
     *  @param values The values to be inserted for the variables
     *  @return The result vector
     */
    template<typename FFTemp>
    std::vector<FFTemp> eval(const std::vector<FFTemp> & values) {
      return static_cast<BlackBoxTemp&>(*this)(values);
    }
    /**
     *  Update internal variables of the black box when the prime field changes.
     *  This function is called from Reconstructor.
     */
    void prime_changed_internal() {
      static_cast<BlackBoxTemp&>(*this).prime_changed();
    }
  };
}