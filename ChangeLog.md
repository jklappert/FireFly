FireFly 2.0.2
=============

New features
------------

 * The factor scan can be performed as a standalone option by calling the
 member function `stop_after_factor_scan` of the `Reconstructor` class.
 The factors can be obtained afterwards with `get_factors_string`.
 An example is given in `example.cpp`.

 * The `ff_insert` executable performs a factor scan only by using the
 `-fs` or `--factorscan` option. The result is written as
 `FUNCTION * FACTOR`.

 * Added on-the-fly and trim functions for the `ShuntingYardParser` class.

 * Added the Meson build system which can be optionally used by the `Kira`
 package. We recommend to build FireFly with `CMake`.

Changes
-------

 * The Horner form generator now rewrites repeated multiplications into an
 exponentiation.

 * When using the option `-m` or `--merge` of the `ff_insert` executable, the
 merged file now has a prefix of the directory that has been used for merging.

 * When using the option `-ni` or `--nointerpolation` of the `ff_insert`
 executable, the `coefficient` directory now carries a suffix of the input file.

 * Improved the memory footprint of the `ShuntingYardParser` class. The RPN is
 only kept in memory on demand when precomputing tokens.

Bug fixes
---------

 * Fixed a bug that prevented the shift from being disabled in rare cases.

 * Add missing operator to the `evaluate` function of the `ShuntingYardParser`
 class. This affacted `^(-` operations only.

 * Fixed a bug in the `ThreadPool` which lead to segmentation faults in
 rare cases.


FireFly 2.0.1
=============

Changes
-------

 * The degree bounds obtained by the shift scan are now used to terminate the
 actual interpolation earlier when possible.

 * Moved the header files to include/firefly.

 * Removed the support for compilable math expressions. The shunting-yard parser
 should be used instead.

Bug fixes
---------

 * Fixed several bugs when restarting from saved states with probes.

 * Improved the CMake script searching for GMP to find installations via the
 distribution. All thanks to Alexander Voigt for solving this issue.

 * Removed trailing whitespaces when building with jemalloc which could lead to
 errors with CMake.

 * The log of the runtime when using the `ff_insert` executable has been fixed.


FireFly 2.0.0
=============

New features
------------

 * Added support for MPI. It can be enabled with the option `-DWITH_MPI=true`.
 The `Reconstructor` class in the master process then requests the `MPIWorker`
 class in the other processes to compute probes. We refer to `README.md` and the
 new `example_mpi.cpp` for more details.

 * Added a scan for univariate factors before the actual interpolation. When
 using the member function `enable_factor_scan` of the `Reconstructor` class,
 FireFly performs the scan and uses the factors for cancellations and
 simplifications.

 * FireFly now compiles an insertion tool for replacement tables in Mathematica
 syntax labeled `ff_insert` by default. This executable serves as an example how
 interpolation algorithms can be applied and can be used to insert integration-
 by-parts tables into amplitudes, for example.

 * FireFly can be linked to the jemalloc library by setting the flag
 `-DWITH_JEMALLOC=true` for CMake. The required paths and libraries are then
 written to firefly.pc and jemalloc is linked to the example and `ff_insert`.

 * The `ShuntingYardParser` can now find same functions. This feature can be
 enabled when constructing a `ShuntingYardParser` with a third argument that is
 set to `true`.

 * The interpolation can now be stopped after the first prime field when the
 member function `reconstruct` of the `Reconstructor` object is called with the
 argument `1`. The intermediate results can be obtained by calling
 `get_result_ff`.

 * A logger file labeled `firefly.log` is now written during the interpolation
 which outputs the current progress and other informational messages.

Changes
-------

 * Changed the `BlackBoxBase` to the curiously recurring template pattern
 (CRTP). This enables FireFly to compute bunches of probes (s. version 1.3.0)
 with a dynamic size depending on the current workload. The bunches are
 implemented in the form of the `FFIntVec` classes, which are static arrays of
 `FFInt` with sizes of powers of 2 up to 128. However, the user has to adopt the
 black box to the CRTP. We refer to `README.me`, `example.cpp`, and the paper of
 version 2.0 for more details.

 * The required C++ standard changed to C++-14.

 * The member function `enable_scan` of the `Reconstructor` class is deprecated.
 It is replaced by the member function `enable_shift_scan`.

 * Added more sophisticated CMake scripts that search for GMP and FLINT.

 * The job scheduling for interpolations in the first prime field or the when
 using the safe mode has been improved.

 * 200 additional prime numbers have been added to the `ReconstHelper` struct.

 * The `Reconstructor` class can now handle empty black boxes.

 * Added pow functions for `FFInt` with a templated power argument. Hence,
 negative integer powers can now be used. We thank Herschel Chawdhry for this
 suggestion.

Bug fixes
---------

 * When generating Horner forms for rational functions with coefficients over Q,
 the coefficient of the maximum degree monomial was incorrect when the numerator
 was 1 and the denominator was not 1. This has been fixed.

 * Fixed `FFInt::pow` yielding a wrong result for powers larger than half of the
 prime number in the default implementation of `FFInt`. We thank Herschel
 Chawdhry for the bug report.


FireFly 1.3.4
=============

Bug fixes
---------

 * When a `PolynomialFF` object contained more than 50 monomials and was altered
 to contain less than 50 monomials, the corresponding Horner form was not
 re-evaluated, thus leading to false evaluations. This has been fixed.


FireFly 1.3.3
=============

Changes
---------

 * Updated PDF.

 * Fixed typos in README.md.

 * Now setting the correct paths to gmp and flint in cmake.

 * Included all necessary flags and paths in firefly.pc.

 * Changed the include path in firefly.pc from ../include/firefly to ../include.


FireFly 1.3.2
=============

Bug fixes
---------

 * When loading from saved states the shift was sometimes disabled only in the
 subsequent prime field leading to avoidable black-box probes. The shift will
 now be disabled in the correct prime field.

 * The parser interpreted some cases of `(+-x^(2*n))` as `+x^(2*n)`. This has
 been fixed. Binary operators were handled correctly.


FireFly 1.3.1
=============

New features
------------

 * When saving the states of reconstruction objects is enabled, every 30 minutes
 all probes which have been calculated are now written to files in the
 `ff_save/probes` directory. Starting from saved states will also read in all
 saved probes and resume the reconstruction from this point. Therefore, one does
 not have to restart from the beginning of the prime field after a crash.

Changes
-------

 * The states of reconstruction objects are now compressed using the `zlib`
 library.

 * Added a new member function for the `Reconstructor` class which reads in a
 whole directory of saved states of reconstruction objects. It can be called
 with `resume_from_saved_state()` when the directory `./ff_save` exist.

 * Added unit tests which check different reconstruction modes. They can be
 executed by calling `make test` in the build directory.

 * Added an implementation of the Xoshiro256** PRNG which is partly used for the
 generation of pseudo random numbers.

 * Minor performance and memory improvements.

Bug fixes
---------

 * Fixed a rare singular-system error which was caused by the sequence of random
 numbers generated by the pcg32 algorithm containing the same number twice. We
 now check the sequences for duplicates and omit them. Thanks to Long Chen for
 providing an example which could be used for finding this error.

 * Fixed a crash with `Nothing left to feed.` which was caused by queuing not
 enough probes in the safe mode with bunch sizes larger than 1.


FireFly 1.3.0
=============

New features
------------

 * Introduced a bunched evaluation of the black box. In addition to the
 `operator()` of `BlackBoxBase` which returns a vector of FFInt, a probe, there
 is now an `operator()` which returns a vector of a vector of FFInt, a vector of
 probes. The default implementation just calls the normal operator several
 times. This feature can help improving the runtime when reaching CPU limits if
 the user can provide an `operator()` which reduces the overhead when it
 computes several probes at once. However, additional threads are always
 preferable to larger bunch sizes.

 * Added a bunched evaluation of parsed functions to the `ShuntingYardParser`.
 It can slightly improve the runtime for large functions.

Changes
-------

 * Changes to the `ShuntingYardParser`:
   - The parsable format changed slightly. The delimiter to mark the end of
   functions `\n` has been replaced by `;`. Spaces and new lines occurring in
   expressions will be removed automatically.

   - The parser now supports unary operators for parentheses. We thank Robert
   Schabinger for this suggestion.

   - The parser now supports negative exponents like `(x+y)^(-10)`. A negative
   exponent has to be used with parentheses.

   - The parser now supports capital letters for variables.

   - The parser now performs a validation of the input by removing white spaces,
   checking parentheses, removing redundant parentheses, and transforming `+-`
   or `-+` to `-`.

 * Without the safe mode, all required probes are now scheduled when changing a
 prime field. This leads to runtime improvements when more than one prime is
 needed.

 * Changed recursive implementation of Thiele and Newton interpolation to
 an iterative one which avoids stack overflows.

 * Added more information to the info messages of the `Reconstructor` class and
 changed them a bit.

 * The safe mode can now handle denominators of monomial coefficients which
 are the prime numbers used for the interpolation.

 * Added a check for the GMP version.

Bug fixes
---------

 * Fixed a rare crash of the `Reconstructor` with `No items to feed anymore`.

 * Fixed another rare crash of the `Reconstructor` with `No items to feed
 anymore`. We thank Mario Prausa for noticing.

 * Fixed 1/0 cases for development versions of FLINT.

 * Fixed the string type for the Intel icpc compiler.


FireFly 1.2.1
=============

Bug fixes
---------

 * Fixed a minor oversampling which occurred in some cases.

 * The `ShuntingYardParser` could parse functions wrongly if the variables
 coincided with internal variables (a,b,c,...,z) which were used for
 replacements. Now, only the user-defined variables are parsed without replacing
 them. We thank Robert Schabinger for noticing.


FireFly 1.2.0
=============

Changes
-------

 * FireFly is now performing a dense-sparse-hybrid approach during the
 interpolation.
 We try to find a sparse shift and then interpolate the rational function densely
 including this shift. When the current highest degree (of either numerator or
 denominator) is fully interpolated, we check if we can use this result to
 interpolate lower degrees sparsely by starting their interpolation again and
 numerically subtracting the effects of the shift if their dense interpolation
 has not finished yet.
 Already (densely or sparsely) interpolated polynomials are removed from the
 homogenized system of equations in `t`.
 Thus, we need the maximal possible number of terms of the rational function
 as black-box probes in the worst case while also utilizing the sparsity of a
 function to some extent. In the average case, this should lead to a drastic
 reduction of the required black-box probes and an overall reduced runtime.

 * When we need a shift in additional prime fields (not in the safe mode), we
 also perform a hybrid approach which first checks which polynomial has the
 maximal number of terms. This polynomial will be calculated sparsely. All
 polynomials with fewer or equal terms than the maximal number (including terms
 generated by a shift) are calculated densely. Polynomials with more terms than
 the maximal number, including terms generated by a shift, are calculated
 sparsely.
 In general, this leads to fewer required black-box probes than a dense approach
 while minimizing the oversampling when a sparse shift is found.

 * `BlackBoxBase` now has a virtual destructor. We thank Philipp Maierhoefer for
 this suggestion.

 * `prime_changed()` in `BlackBoxBase` now has a default implementation which
 does nothing instead of requiring the user to define it. We thank Mario Prausa
 for this suggestion.

 * `enable_scan()` in the safe mode now just throws a warning instead of exiting
 the program. We thank Mario Prausa for this suggestion.

Bug fixes
---------

 * The `Reconstructor` class now has a destructor which deletes all dynamically
 allocated memory.

 * The safe mode is usable again. We thank Mario Prausa for noticing the
 problems.

 * If the number of variables was set to 1 and the black box was just a constant,
 the interpolation failed. This has been fixed and constants for one variable
 can be interpolated again. We thank Mario Prausa for noticing.

 * The `ShuntingYardParser` could not parse negative variables. This has been fixed.

 * The `ShuntingYardParser` could not handle global signs. This has been fixed.

 * Fixed a deadlock when starting from a saved state.

 * `set_tags()` and `resume_from_saved_state()` now work more consistently.

 * `-lflint` is only added to the pkgconfig file when FLINT is actually used.


FireFly 1.1.1
=============

Bug fixes
---------

 * The `Reconstructor` class can reconstruct zeros again.


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

 * Added a "safe mode" which interpolates the black box from scratch over
 each prime field to be sensitive to unlucky primes, zeros, and other
 errors. It can be used by calling the member function `set_safe_interpolation()`
 of the `Reconstructor` class.

 * The black box is no longer a function of `Reconstructor` but instead a functor.
 `source/include/Reconstructor.hpp` contains the base class `BlackBoxBase`.
 The user has to define his black box as a derived class and provide a constructor,
 the evaluation `std::vector<FFInt> operator()(const std::vector<FFInt>& values)`,
 and the function `void prime_changed()` which allows the user to change class
 variables when the prime field changes. `Reconstructor` has to be initialized with
 a reference to the black box. The example class `BlackBoxUser` can be found in
 `example.cpp`

 * Added a script which converts a list of functions in Mathematica syntax
 to compilable C++ code to perform interpolations with them. This can be
 helpful for arithmetic with functions with many terms where other programs
 fail. The script can be found in the `mma_2_ff` directory.

 * Added a shunting-yard parser to parse a collection of rational functions
 for functional evaluation. This skips the compilation steps of the Mathematica
 to C++ conversion script.

 * Added a Horner representation of rational functions which may result in
 faster evaluation times when needed. Check out `source/include/HornerGenerator.hpp`
 and the new member functions of `RationalFunctions` and `Polynomial`, respectively.

 * Added dense algorithms for matrix manipulation, e.g., computing the inverse,
 the determinant, or the solution of a system of equations. Additionally,
 LU decompositions are supported. Further information can be found in
 `source/include/DenseSolver.hpp`.

Changes
-------

 * Generate 64-bit anchor points and shifts instead of 32-bit.

 * Many small runtime improvements and minor bug fixes.

 * Discard prime fields over which the black box evaluates to zero in the first
 probe. This makes the code safer regarding unlucky primes and zeros.

Bug fixes
---------

 * When the member function `to_string` of the `RationalFunction` class was called
 on a zero polynomial an error was thrown which was actually not an error. Thanks
 to Long Chen.


FireFly 1.0.0
=============

This is the initial release of the FireFly library. FireFly interpolates rational
functions over finite (prime) fields Z_p.
