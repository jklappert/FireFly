#pragma once

#include "FFInt.hpp"

namespace firefly {

  typedef std::vector<std::vector<FFInt>> mat_ff;
  /**
   *  Calculates the inverse of a matrix using Gauss-Jordan
   * 
   *  @param a input matrix build of FFInts
   *  @param n the size of a
   */
  void calc_inverse(mat_ff& a, uint32_t n);
  /**
  *  Solves the given system of equations using a Gauss-Jordan algorithm
  *  @param a the (n x (n + 1)) matrix which represents the system of equations, e.q., (x, x^2, | f(x))
  *  @param n the number of equations
  *  @return The solved coefficients of the systems
  */
  std::vector<FFInt> solve_gauss_system(mat_ff& a, uint32_t n);
  /**
   *  Solves a system of equations, A*x=b, using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param b the right hand side of the system
   *  @param n the size of a
   *  @return the result vector x
   */
  std::vector<FFInt> solve_lu(mat_ff& a, std::vector<int>& p, const std::vector<FFInt>& b, uint32_t n);
  /**
   *  Inverts a matrix using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param ia will be filled by the inverse of a
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param n the size of a
   */
  void calc_inverse_lu(const mat_ff& a, mat_ff& ia, std::vector<int>& p, uint32_t n);
  /**
   *  Calculates the determinant of a matrix using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param n the size of a
   *  @return the determinant of a
   */
  FFInt calc_determinant_lu(const mat_ff& a, std::vector<int>& p, uint32_t n);
  /**
   *  Decomposes a matrix accodring to LU decomposition and saves its form in the given input
   *  @param a the matrix of which a LU decomposition should be performed. The result is saved in a.
   *  @param p the permutation matrix stored as an vector of size n+1 containing column indices, where the permutation matrix has 1.
   *  @param n the size of the matrix
   */
  void calc_lu_decomposition(mat_ff& a, std::vector<int>& p, uint32_t n);

}
