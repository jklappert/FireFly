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

//#include "FFInt.hpp" // for examples
#include "Logger.hpp"

#include <vector>

namespace firefly {

  // Some examples of the dense solver functions

  // Initialize matrices
  /*mat_ff<FFInt> a = {{0*33,17,25},{59595, 989983749,99},{23213213, 4354354353,0*434232324}};
  mat_ff<FFInt> a_2 = a;
  mat_ff<FFInt> a_3 = a;
  a = {{0*33,17,25,10},{59595, 989983749,99,14},{23213213, 4354354353,0*434232324,200}};
  mat_ff<FFInt> inv;
  //----------------------------------------------------------------------------
  // Calculate inverse of a using Gauss-Jordan algorithm
  calc_inverse(a_3, 3);
  // Print inverse of a_3 obtained with Gauss-Jordan
  std::cout << "Inverse\n" << a_3[0][0] << " " << a_3[0][1] << " " << a_3[0][2] << "\n"
  << a_3[1][0] << " " << a_3[1][1] << " " << a_3[1][2] << "\n"
  << a_3[2][0] << " " << a_3[2][2] << " " << a_3[2][2] << "\n";
  // Solve a using Gauss algorithm
  std::vector<FFInt> res = solve_gauss_system(a, 3);
  std::cout << "Sol Gauss\n" << res[0] << " " << res[1] << " " << res[2] << "\n";
  //----------------------------------------------------------------------------
  // Initialize permutation matrix for LU decompositions
  std::vector<int> p;
  // Initialize vector for solution of LU system
  std::vector<FFInt> b = {10,14,200};
  // Calculate LU decomposition of a_2
  calc_lu_decomposition(a_2, p, 3);
  // Solve a_2 using LU decomposition
  res = solve_lu(a_2, p, b, 3);
  std::cout << "Sol LU\n" << res[0] << " " << res[1] << " " << res[2] << "\n";
  // Calculate inverse of a_2 using LU decomposition
  calc_inverse_lu(a_2, inv, p, 3);
  // Calculate determinant of a_2 using LU decomposition
  std::cout << "Det LU " << calc_determinant_lu(a_2, p, 3) << "\n";
  // Print inverse obtained using LU decomposition
  std::cout << "Inverse LU\n" << inv[0][0] << " " << inv[0][1] << " " << inv[0][2] << "\n"
  << inv[1][0] << " " << inv[1][1] << " " << inv[1][2] << "\n"
  << inv[2][0] << " " << inv[2][1] << " " << inv[2][2] << "\n";*/
//TODO modify result above
  /*      //std::vector<FFIntTemp> b(1);
      //solve_lu(mat, p, b, 2);

      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //calc_inverse(mat, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //solve_gauss_system(mat, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //p = {};
      //mat_ff<FFIntTemp> mat2;
      //calc_inverse_lu(mat, mat2, p, 2);
      //mat = {{result[0], result[1]}, {result[2], result[3]}};
      //p = {};
      //calc_determinant_lu(mat, p, 2);*/

  template<typename FFIntTemp>
  using mat_ff = std::vector<std::vector<FFIntTemp>>;

  /**
   *  Calculates the inverse of a matrix using Gauss-Jordan
   *  @param a input matrix build of FFInts
   *  @param n_ the size of a
   */
  template<typename FFIntTemp>
  void calc_inverse(mat_ff<FFIntTemp>& a, uint32_t n_) {
    int n = static_cast<int>(n_);

    // Augment a with unit matrix
    for (int i = 0; i < n; ++i) {
      std::vector<FFIntTemp> dum(n);
      dum[i] = FFIntTemp(1);
      std::vector<FFIntTemp> tmp = a[i];
      tmp.insert(tmp.end(), dum.begin(), dum.end());
      a[i] = tmp;
    }

    for (int i = 0; i < n; ++i) {
      // search for maximum in this column
      FFIntTemp max_el = a[i][i];
      int max_row = i;

      for (int k = i + 1; k < n; ++k) {
        if (a[k][i] > max_el) {
          max_el = a[k][i];
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (int k = i; k < 2 * n; ++k) {
        FFIntTemp tmp = a[max_row][k];
        a[max_row][k] = a[i][k];
        a[i][k] = tmp;
      }

      // make all rows below this one 0 in current column
      for (int k = i + 1; k < n; ++k) {
        FFIntTemp c = -a[k][i] / a[i][i];

        for (int j = i; j < 2 * n; ++j) {
          if (i == j)
            a[k][j] = FFIntTemp(0);
          else
            a[k][j] += c * a[i][j];
        }
      }
    }

    // solve equation ax=b for an upper triangular matrix a
    mat_ff<FFIntTemp> res(n, std::vector<FFIntTemp>(n));

    for (int i = n - 1; i >= 0; i--) {
      for (int k = n; k < 2 * n; ++k) {
        a[i][k] /= a[i][i];
      }

      for (int row_mod = i - 1; row_mod >= 0; row_mod--) {
        for (int column_mod = n; column_mod < 2 * n; ++column_mod) {
          a[row_mod][column_mod] -= a[i][column_mod] * a[row_mod][i];
        }
      }

      // Remove the unit matrix from the result
      for (int k = n; k < 2 * n; ++k) {
        res[i][k - n] = a[i][k];
      }
    }

    a = res;
  }
  /**
  *  Solves the given system of equations using a Gauss-Jordan algorithm
  *  @param a the (n x (n + 1)) matrix which represents the system of equations, e.q., (x, x^2, | f(x))
  *  @param n the number of equations
  *  @param force forces to give solution also for singular systems
  *  @return The solved coefficients of the systems
  */
 template<typename FFIntTemp>
  std::vector<FFIntTemp> solve_gauss_system(mat_ff<FFIntTemp>& a, uint32_t n, bool force = false) {
    // Transform the matrix in upper triangular form
    for (uint32_t i = 0; i < n; ++i) {
      // search for maximum in this column
      FFIntTemp max_el = a[i][i];
      uint32_t max_row = i;

      for (uint32_t k = i + 1; k < n; ++k) {
        FFIntTemp tmp = a[k][i];

        if (tmp > max_el) {
          max_el = tmp;
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (uint32_t k = i; k < n + 1; ++k) {
        FFIntTemp tmp = a[max_row][k];
        a[max_row][k] = a[i][k];
        a[i][k] = tmp;
      }

      // Make all rows below this one zero in the current column
      for (uint32_t k = i + 1; k < n; ++k) {
        FFIntTemp c = -a[k][i] / a[i][i];

        for (uint32_t j = i; j < n + 1; ++j) {
          if (i == j) a[k][j] = FFIntTemp(0);
          else a[k][j] += c * a[i][j];
        }
      }
    }

    std::vector<FFIntTemp> results(n);

    if (force || a[n - 1][n - 1] != FFIntTemp(0)) { // TODO new == operator?
      // Solve equation A * x = b for an upper triangular matrix
      for (int i = n - 1; i >= 0; i--) {
        results[i] = a[i][n] / a[i][i];

        for (int k = i - 1; k >= 0; k--) {
          a[k][n] -= a[k][i] * results[i];
        }
      }
    } else {
      /*for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
          std::cout << a[i][j] << " ";
        }

        std::cout << "\n";
      }*/

      ERROR_MSG("Singular system of equations!");
      std::exit(EXIT_FAILURE);
    }

    return results;
  }
  /**
   *  Solves a system of equations, A*x=b, using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param b the right hand side of the system
   *  @param n the size of a
   *  @return the result vector x
   */
  template<typename FFIntTemp>
  std::vector<FFIntTemp> solve_lu(mat_ff<FFIntTemp>& a, std::vector<int>& p, const std::vector<FFIntTemp>& b, uint32_t n) {
    std::vector<FFIntTemp> x(n);

    for (uint32_t i = 0; i < n; ++i) {
      x[i] = b[p[i]];

      for (uint32_t k = 0; k < i; ++k)
        x[i] -= a[i][k] * x[k];
    }

    for (int i = n - 1; i >= 0; i--) {
      for (int k = i + 1; k < static_cast<int>(n); ++k)
        x[i] -= a[i][k] * x[k];

      x[i] = x[i] / a[i][i];
    }

    return x;
  }
  /**
   *  Inverts a matrix using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param ia will be filled by the inverse of a
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param n the size of a
   */
  template<typename FFIntTemp>
  void calc_inverse_lu(const mat_ff<FFIntTemp>& a, mat_ff<FFIntTemp>& ia, const std::vector<int>& p, uint32_t n) {
    ia = mat_ff<FFIntTemp>(n, std::vector<FFIntTemp> (n));

    for (uint32_t j = 0; j < n; ++j) {
      for (uint32_t i = 0; i < n; ++i) {
        if (p[i] == static_cast<int>(j))
          ia[i][j] = FFIntTemp(1);
        else
          ia[i][j] = FFIntTemp(0)
          ;

        for (uint32_t k = 0; k < i; ++k)
          ia[i][j] -= a[i][k] * ia[k][j];
      }

      for (int i = n - 1; i >= 0; i--) {
        for (int k = i + 1; k < static_cast<int>(n); ++k)
          ia[i][j] -= a[i][k] * ia[k][j];

        ia[i][j] = ia[i][j] / a[i][i];
      }
    }
  }
  /**
   *  Calculates the determinant of a matrix using LU factorization
   *  @param a input matrix build of FFInts and is already LU decomposed
   *  @param p the permutation matrix obtained during the LU decomposition of a
   *  @param n the size of a
   *  @return the determinant of a
   */
  template<typename FFIntTemp>
  FFIntTemp calc_determinant_lu(const mat_ff<FFIntTemp>& a, const std::vector<int>& p, uint32_t n) {
    FFIntTemp det = a[0][0];

    for (uint32_t i = 1; i < n; ++i)
      det *= a[i][i];

    if ((p[n] - n) % 2 == 0)
      return det;
    else
      return -det;
  }
  /**
   *  Decomposes a matrix accodring to LU decomposition and saves its form in the given input
   *  @param a the matrix of which a LU decomposition should be performed. The result is saved in a.
   *  @param p the permutation matrix stored as an vector of size n+1 containing column indices, where the permutation matrix has 1.
   *  @param n the size of the matrix
   */
  template<typename FFIntTemp>
  void calc_lu_decomposition(mat_ff<FFIntTemp>& a, std::vector<int>& p, uint32_t n) {
    uint32_t i, k, j, max_row;
    FFIntTemp max_el;
    std::vector<FFIntTemp> tmp;
    p = std::vector<int> (n + 1);

    for (i = 0; i <= n; ++i)
      p[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < n; ++i) {
      max_el = FFIntTemp(0);
      max_row = i;

      for (k = i; k < n; ++k)
        if (a[k][i] > max_el) {
          max_el = a[k][i];
          max_row = k;
        }

      if (max_row != i) {
        //pivoting p
        j = p[i];
        p[i] = p[max_row];
        p[max_row] = j;

        //pivoting rows of A
        tmp = a[i];
        a[i] = a[max_row];
        a[max_row] = tmp;

        //counting pivots starting from N (for determinant)
        p[n]++;
      }

      for (j = i + 1; j < n; ++j) {
        if (a[i][i] == FFIntTemp(0)) {// TODO new == operator?
          ERROR_MSG("Division by zero while calculating LU decomposition.");
          std::runtime_error("Division by zero while calculating LU decomposition.");
        }

        a[j][i] /= a[i][i];

        for (k = i + 1; k < n; ++k)
          a[j][k] -= a[j][i] * a[i][k];
      }
    }
  }

}
