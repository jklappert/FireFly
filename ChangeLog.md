FireFly 1.1.0
=============

New features
------------
 * Added Ben-Or/Tiwari sparse univariate polynomial interpolation algorithm
 and combined it with the dense Newton interpolation algorithm to a
 racing algorithm during Zippel's algorithm. This reduces the number of
 black-box probes for sparse multivariate polynomials significantly.
 Thanks to Sven Yannick Klein.

 * During the interpolation of the univariate rational auxiliary functions
 constants in numerator or denominator which are not used for normalization
 are removed from the systems of equations to reduce the number of black-box
 probes.

 * Added a "save mode" which interpolates the black box from scratch over
 each prime field to be sensitive to unlicky primes, zeros and other
 errors. It can be used by calling the member function `set_save_interpolation()`
 of the `Reconstructor` class.

 * Added a script which converts a list of functions in Mathematica syntax
 to compilable C++ code to perform interpolations with them. This can be
 helpful for arithmetic with functions with many terms where other programs
 fail. The script can be found in the `mma_2_ff` directory.

 * Added dense algorithms for matrix manipulation, e.g., computing the inverse,
 the determinant, or solutions for system of equations. Additionaly,
 LU decompositions are supported. Further information can be found in
 `source/DenseSolver.hpp`.

Changes
-------

 * Generate 64-bit anchor points and shifts instead of 32-bit.

 * Many small runtime improvements and minor bug fixes.


FireFly 1.0.0
=============

This is the initial release of the FireFly library. The FireFly
interpolates rational functions over finite fields Z_p.
