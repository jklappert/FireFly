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

#pragma once

//#include "FFInt.hpp" // for examples
#include "firefly/Logger.hpp"

#include <vector>

namespace firefly {
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
