#pragma once

#include "FFInt.hpp"

namespace firefly {

  typedef std::vector<std::vector<FFInt>> mat_ff;
  /**
   *  Calculates the cofactor of a[p][q] of tmp[][]
   *  @param a input matrix build of FFInts
   *  @param tmp a temporary object needed for the cofactors. It's a matrix build from FFInts
   *  @param p current row
   *  @param q current column
   *  @param n current rank of a
   */
  void calc_cofactor(const mat_ff& a, mat_ff& tmp, uint32_t p, uint32_t q, uint32_t n);

  /**
   *  Calculates the determinat of a matrix
   *  @param a input matrix build of FFInts
   *  @param n rank of a
   */
  FFInt calc_determinant(const mat_ff& a, uint32_t n);

  /**
   *  Calculates the adjoint of a matrix
   *  @param a input matrix build of FFInts
   *  @param adj the adjoint matrix calculated by this function
   *  @param n the size of a
   */
  void calc_adjoint(const mat_ff& a, mat_ff& adj, uint32_t n);

  /**
   *  Calculates the inverse of a matrix
   *  @param a input matrix build of FFInts
   *  @param inv inverse matrix to a if a is invertible
   *  @param n the size of a
   *  @return True if a is invertible
   */
  bool calc_inverse(const mat_ff& a, mat_ff& inv, uint32_t n);
  
  /**
   *  Calculates the inverse of a matrix using Gauss-Jordan
   * 
   *  @param a input matrix build of FFInts
   *  @param n the size of a
   */
  void calc_inverse_2(mat_ff& a, uint32_t n);

  /**
  *  Solves the given system of equations using a Gauss-Jordan algorithm
  *  @param coef_mat the matrix which represents the system of equations
  *  @param num_eqn the number of equations
  *  @return The solved coefficients of the systems
  */
  std::vector<FFInt> solve_gauss_system(mat_ff& coef_mat, uint32_t num_eqn);
}
