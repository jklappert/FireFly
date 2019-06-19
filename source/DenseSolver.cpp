#include "DenseSolver.hpp"
#include "Logger.hpp"

namespace firefly {

  void calc_cofactor(const mat_ff& a, mat_ff& tmp, uint32_t p, uint32_t q, uint32_t n) {
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; ++row) {
      for (int col = 0; col < n; ++col) {
        //  Copying into temporary matrix only those element
        //  which are not in given row and column
        if (row != p && col != q) {
          tmp[i][j++] = a[row][col];

          // Row is filled, so increase row index and
          // reset col index
          if (j == n - 1) {
            j = 0;
            ++i;
          }
        }
      }
    }
  }

  FFInt calc_determinant(const mat_ff& a, uint32_t n) {
    FFInt det = 0;

    //  Base case : if matrix contains single element
    if (n == 1)
      return a[0][0];

    mat_ff tmp (n, std::vector<FFInt> (n)); // To store cofactors

    int sign = 1;  // To store sign multiplier

    // Iterate for each element of first row
    for (uint32_t i = 0; i < n; ++i) {
      // Getting cofactor of a[0][i]
      calc_cofactor(a, tmp, 0, i, n);
      det += sign * a[0][i] * calc_determinant(tmp, n - 1);

      // terms are to be added with alternate sign
      sign = -sign;
    }

    return det;
  }

  void calc_adjoint(const mat_ff& a, mat_ff& adj, uint32_t n) {
    if (n == 1) {
      adj[0][0] = 1;
      return;
    }

    // temp is used to store cofactors of a[][]
    int sign = 1;
    mat_ff tmp (n, std::vector<FFInt> (n));

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        // Get cofactor of a[i][j]
        calc_cofactor(a, tmp, i, j, n);

        // sign of adj[j][i] positive if sum of row
        // and column indexes is even.
        sign = ((i + j) % 2 == 0) ? 1 : -1;

        // Interchanging rows and columns to get the
        // transpose of the cofactor matrix
        adj[j][i] = (sign) * (calc_determinant(tmp, n - 1));
      }
    }
  }

  bool calc_inverse(const mat_ff& a, mat_ff& inv, uint32_t n) {
    // Find determinant of a[][]
    FFInt det = calc_determinant(a, n);

    if (det == 0) {
      ERROR_MSG("Singular matrix, cannot find its inverse.");
      return false;
    }

    // Find adjoint
    mat_ff adj (n, std::vector<FFInt> (n));

    calc_adjoint(a, adj, n);

    // Clear inverse
    inv = mat_ff (n, std::vector<FFInt> (n));

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        inv[i][j] = adj[i][j] / det;

    return true;
  }

  std::vector<FFInt> solve_gauss_system(mat_ff& coef_mat, uint32_t num_eqn) {
    // Transform the matrix in upper triangular form
    for (uint32_t i = 0; i < num_eqn; ++i) {
      // search for maximum in this column
      FFInt max_el = coef_mat[i][i];
      uint32_t max_row = i;

      for (uint32_t k = i + 1; k < num_eqn; ++k) {
        FFInt tmp = coef_mat[k][i];

        if (tmp.n > max_el.n) {
          max_el = tmp;
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (uint32_t k = i; k < num_eqn + 1; ++k) {
        FFInt tmp = coef_mat[max_row][k];
        coef_mat[max_row][k] = coef_mat[i][k];
        coef_mat[i][k] = tmp;
      }

      // Make all rows below this one zero in the current column
      for (uint32_t k = i + 1; k < num_eqn; ++k) {
        FFInt c = -coef_mat[k][i] / coef_mat[i][i];

        for (uint32_t j = i; j < num_eqn + 1; ++j) {
          if (i == j) coef_mat[k][j] = FFInt(0);
          else coef_mat[k][j] += c * coef_mat[i][j];
        }
      }
    }

    std::vector<FFInt> results(num_eqn);

    if (coef_mat[num_eqn - 1][num_eqn - 1] != 0) {
      // Solve equation A * x = b for an upper triangular matrix
      for (int i = num_eqn - 1; i >= 0; i--) {
        results[i] = coef_mat[i][num_eqn] / coef_mat[i][i];

        for (int k = i - 1; k >= 0; k--) {
          coef_mat[k][num_eqn] -= coef_mat[k][i] * results[i];
        }
      }
    } else {
      ERROR_MSG("Singular system of equations!");
      throw std::runtime_error("Singular system of equations!");
    }

    return results;
  }

}
