# FireFly

[![Build Status](https://travis-ci.org/jklappert/FireFly.svg?branch=master)](https://travis-ci.org/jklappert/FireFly) [![](https://img.shields.io/github/tag/jklappert/firefly)](https://gitlab.com/firefly-library/firefly/-/tags/2.0.2)

FireFly is a reconstruction library for rational functions written in C++.

Please refer to these papers when using FireFly:
* [1] J. Klappert and F. Lange, *Reconstructing Rational Functions with FireFly*, [[*Comput.Phys.Commun.* **247** (2020) 106951](https://doi.org/10.1016/j.cpc.2019.106951)], [[1904.00009](https://arxiv.org/abs/1904.00009)]
* [2] J. Klappert, S.Y. Klein, and F. Lange, *Interpolation of Dense and Sparse Rational Functions and other Improvements in FireFly*, [[2004.01463](https://arxiv.org/abs/2004.01463)]

## Table of contents
* [Requirements](#requirements)
* [Building FireFly](#building-firefly)
* [Reconstructing functions](#reconstructing-functions)
* [Directly parse collections of rational functions](#directly-parse-collections-of-rational-functions)
* [Code Documentation](#code-documentation)
* [Compiling Code with FireFly](#compiling-code-with-firefly)
* [FireFly Executable (in Development)](#firefly-executable-in-development)

## Requirements
FireFly requires:
* C++ compiler supporting C++14
* [CMake](https://cmake.org/) >= 3.1
* [FLINT](http://www.flintlib.org/) >= 2.5 (optional)
* [GMP](https://gmplib.org/) >= 6.1.2 (compiled with --enable-cxx)
* [zlib](https://www.zlib.net/) >= 1.2.11
* MPI >= 3 (optional, mpich recommended)

#####  Third party code
FireFly uses the following third party code:
* [tinydir](https://github.com/cxong/tinydir)
* [gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/)

#####  Implemented interpolation algorithms
FireFly implements the following third party interpolation algorithms:
* Ben-Or/Tiwari algorithm [[10.1145/62212.62241](https://dx.doi.org/10.1145/62212.62241)]
* Newton interpolation algorithm
* Racing algorithm [[10.1145/62212.62241](https://dx.doi.org/10.1145/345542.345629)], [[10.1016/S0747-7171(03)00088-9](https://dx.doi.org/10.1016/S0747-7171(03)00088-9)]
* Sparse rational function interpolation algorithm [[10.1016/j.tcs.2010.11.050](https://dx.doi.org/10.1016/j.tcs.2010.11.050)]
* Thiele interpolation algorithm
* Zippel algorithm [[10.1016/S0747-7171(08)80018-1](https://dx.doi.org/10.1016/S0747-7171(08)80018-1)]

## Building FireFly
FireFly uses CMake to generate files for build automation. To build FireFly one should first create a separate `build` directory inside FireFly's top directory. Afterwards, `cmake` should be called:

```
cd $FIREFLY_PATH
mkdir build
cd build
cmake -DWITH_FLINT=true .. # Without FLINT: -DWITH_FLINT=false or omit this
```

After calling `cmake` the build directory contains all required build files. Assuming that GNU make is used, one can start the build by running

```
make
```

and afterwards install it with

```
make install
```

By default the code is compiled with optimizations.

If FLINT is used for modular arithmetic and it cannot be found in the default system directories, one has to add the additional flags:

```
-DFLINT_INCLUDE_DIR=$FLINT_INC_PATH -DFLINT_LIBRARY=$FLINT_LIB_PATH
```

where `FLINT_LIB_PATH` is the absolute path pointing to the shared library of FLINT.

If the GMP version installed in the system directories does not match 6.1.2 or CMake does not find GMP, the paths for GMP can be set with the flags:
```
-DGMP_INCLUDE_DIRS=$GMP_INC_PATH -DGMP_LIBRARIES=$GMP_LIB_PATH
```
where `GMP_INC_PATH` is the absolute path to the directory of where the include files can be found (`gmpxx.h` is required) and `GMP_LIB_PATH` is the absolute path to the GMP library.

One can include the support of the `jemalloc` memory allocation library by adding the flag
```
-DWITH_JEMALLOC=true
```
after the `jemalloc` library is installed and its executables can be found.

The MPI version of FireFly can be built by adding the flag
```
-DWITH_MPI=true
```
If the MPI version found with CMake is not satisfactory, one can set custom paths to the include directory and library as:
```
-DMPI_CXX_INCLUDE_PATH=$MPI_INC_PATH -DMPI_CXX_LIBRARIES=$MPI_LIB_PATH
```
where `MPI_INC_PATH` is the absolute path to the directory where the include files can be found and `MPI_LIB_PATH` is the absolute path to the shared library (`libmpi.so`, `libmpich.so`, ...) of your favored MPI implementation.

Below is a list of all build options for `FireFly`:

* `-DWITH_FLINT=true` (default: `false`): Employ FLINT for the modular arithmetic (highly recommended).
* `-DCUSTOM=true` (default: `false`): Employ a custom modular arithmetic.
* `-DWITH_JEMALLOC=true` (default: `false`): Employ jemalloc for malloc operations (highly recommended).
* `-DWITH_MPI=true` (default: `false`): Enable MPI support.
* `-DGMP_INCLUDE_DIRS=$GMP_INC_PATH` (default: not set): If GMP is not automatically found, one has to manually set the absolute path to its header files (`gmpxx.h` is required).
* `-DGMP_LIBRARIES=$GMP_LIB_PATH` (default: not set): If GMP is not automatically found, one has to manually set the absolute path to the shared library.
* `-DFLINT_INCLUDE_DIR=$FLINT_INC_PATH` (default: not set): If FLINT is not automatically found, one has to manually set the absolute path to its header files.
* `-DFLINT_LIBRARY=$FLINT_LIB_PATH` (default: not set): If FLINT is not automatically found, one has to manually set the absolute path to the shared library.
* `-DMPI_CXX_INCLUDE_PATH=$MPI_INC_PATH` (default: not set): If MPI is not automatically found, one has to manually set the absolute path to its header files.
* `-DMPI_CXX_LIBRARIES=$MPI_LIB_PATH` (default: not set): If MPI is not automatically found, one has to manually set the absolute path to the shared library.
* `-DBUILD_BENCH=true` (default: `false`): Build the benchmark executables `benchmarks` and `benchmarks_no_bt`, which perform some of the benchmarks from Ref. [2].
* `-DFIREFLY_EXECUTABLE=true` (default: `false`): Build and install the executable `firefly`, cf. [FireFly Executable (in Development)](#firefly-executable-in-development).

Important standard options of CMake are:

* `-DCMAKE_INSTALL_PREFIX=$INSTALL_PATH` (default: standard path depending on the system): Path to which FireFly is installed.
* `-DCMAKE_BUILD_TYPE=Debug` (default: `Release`): Compile with warnings in debug mode.

## Reconstructing functions
To reconstruct functions with FireFly, it offers an interface which directly makes use of a thread pool for the parallel reconstruction of various functions over the same prime field. Additionally, black-box probes are calculated in parallel.

The black box is implemented as functor following the curiously recurring template pattern. The user has to define the black box as a derived class of `BlackBoxBase` and provide a constructor and the evaluation of the black box:

```cpp
class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
public:
  BlackBoxUser(...) {
    ...
  }
  template<typename FFIntTemp>
  std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values) {
    ...
  }
}
```

Optionally, the user can provide the function `void prime_changed()` which allows the user to change member variables when the prime field changes.

The actual reconstruction is done by the `Reconstructor` object which is initialized with the number of variables `n_var`, the number of threads `n_thr` that should be used, and the `BlackBoxUser` object `bb`:

```cpp
BlackBoxUser bb(...);
Reconstructor<BlackBoxUser> rec(n_var, n_thr, bb);
```

The reconstruction is then started by

```cpp
rec.reconstruct();
```

The reconstruction will run from this point until it is finished. The results can be obtained with

```cpp
rec.get_result();
```

We refer to the `example.cpp` file and the code documentation for additional options and features.


## Directly parse collections of rational functions
FireFly provides a parser class for rational functions. The functions will be stored in reverse polish notation to be evaluated for a given parameter point. Parsing collections of rational functions can be done with the `ShuntingYardParser` class. It has to be constructed with a path to a file which contains the rational functions and a vector of the occurring variables:

```cpp
ShuntingYardParser parser(path, vars);
```

Here, `path` is a string containing the path to the file in which the functions are stored and `vars` is a vector of strings of the variables. The collection of functions have to be separated by a semicolon `;` to be identified correctly, e.g.

```
(x+y*3)/(z^2+x);
(12132132323213213212/33*x + 12 * 3 - x^100*y^2)/(3*y^5);
```

The corresponding vector `vars` would thus be

```
std::vector<std::string> vars = {"x","y","z"};
```

Variables are limited to a length of at most 16 characters and can consist of lower and upper case letters, i.e. `a,...,z` and `A,...,Z`, and numbers, e.g. `s12`.

The functions have to be parsed only once and can be evaluated afterwards calling

```cpp
parser.evaluate(values);
// Evaluates the black-box functions with precomputed values which is faster than evaluate().
// Its usage requires parser.precompute_tokens() after the field has changed.
//parser.evaluate_pre(values);
```

where `values` is a vector which contains the parameter point at which the functions should be evaluated. The function `evaluate` returns a vector of `FFIntTemp` objects, i.e. `FFInt` or `FFIntVec`, which is filled with the values of the evaluated functions in the same order as the functions are defined in the input file. Thus, it can be directly used in the `BlackBox` functor of FireFly. An example file is given in `parser_test/s_y_4_v.m`. Note that only the operators

```
+, -, *, /, ^
```

are supported. Negative exponents like `x^(-10)` have to be set in parentheses. For improved runtime, the evaluation should be performed with `evaluate_pre` which uses precomputed values of the monomial coefficients for the current field.

For convenience, FireFly also provides a script which converts a list of rational functions (stored as an expression list of Mathematica) to FireFly's parsable format. It is located in the `mma_2_ff` directory and can be executed with

```
./convert_to_sy.sh $FILE
```

where `$FILE` contains the list of rational functions.

**The following operations are not supported:**

* Any kind of implicit operators like `3 x`. This should read `3*x` instead.
* Negative exponents without parentheses like `x^-5`. This should read `x^(-5)` instead.
* Unevaluated exponents like `x^(3+7)`. This should read `x^10` instead.
* Operators followed by operators should be separated, i.e., for example, `3*-x` should read `3*(-x)` or `-3*x`. Only `+-` or `-+` will be interpreted as `-`.

## Code Documentation
Doxygen can be used to generate code documentation. To generate the documentation, run

```
make doc
```

The generated documentation can be found in `doc/html/index.html` or in `doc/latex`. Using the latter directory, calling `make` will create a pdf file.

## Compiling Code with FireFly

In order to successfully compile and link another program with FireFly, one has to set the correct compiler and linker flags.
These can be found in the file `firefly.pc`, which is installed together with FireFly into the subdirectory `lib/pkgconfig` of the path set by the user with `-DCMAKE_INSTALL_PREFIX` (or the default path on the system).
The easiest possibility to compile and link code with FireFly is to import it employing the `pkg-config` support of a build system.
When you do not use a build system, the required flags can be extracted from `firefly.pc` using `pkg-config`, e.g.

```
FF_CFLAGS=$(pkg-config --cflags firefly)
FF_LIBS=$(pkg-config --libs firefly)
```

Then, the code can be compiled and linked with

```
c++ $FF_CFLAGS <USER_CODE> $FF_LIBS
```

If `pkg-config` is not installed on the system, the flags can simply be read off `firefly.pc`.

## FireFly Executable (in Development)

`FireFly` offers the executable `firefly`, which can read in probes from a file and interpolate and reconstruct the functions from it.
This feature is still in development and not officially released.
Currently, only the most basic options for `FireFly` are supported, i.e. neither the scan for a sparse shift nor the scan for univariate factors are performed.

With the CMake option `-DFIREFLY_EXECUTABLE=true`, one can build and install the executable if neither a custom modular arithmetic nor `MPI` is used.

The user has to set the anchor points and the shifts for all variables in all prime fields and write them to the files `anchor_points` and `shifts` in the directory, in which `firefly` is supposed to run.
The anchor points and shifts are related to the values for the variables through

```
z_{i, j_i, k} = t_k * a_i^{j_i} + s_i
```

where `t_k` is an arbitrary random number, `a_i` is the anchor point for variable `i`, `s_i` is the shift for variable `i`, and `j_i > 0` is some power for the anchor point `a_i`.
In the publications on `FireFly`, the powers are referred to as `zi_orders`.
The anchor points for the first variable are always set to 1 and they do not have to be set by the user.
FireFly requests each probe to be computed at some tuple of the powers `j_i`, e.g. the tuple `(3,1)` translates to

```
z_{1, x, k} = t_k * 1^x + s_1 , z_{2, 3, k} = t_k * a_2^3 + s_2 , z_{3, 1, k} = t_k * a_3^1 + s_3 ,
```

where `x` is irrelevant, because the anchor point is 1.

This special form is required for two algorithms used in polynomial interpolation, namely solving Vandermonde systems (see the last two pages of Sct. 2.1 in Ref. [1]) and the Ben-Or/Tiwari algorithm (see Sct. 2.1.1 in Ref. [2]).
Usually, the powers `j_i` monotonically increase in steps of 1 until the algorithms terminate with the current variable, e.g. a possible order would be `(1,1), (2,1), ..., (10,1), (1,2), (2,2), (3,2), ...`

For three variables, the file `anchor_points` may then look like

```
6747833927250002263 1719011221300048416
6706258047612828567 5195796091462229223
```

The anchor points may coincide between different prime fields, but we encourage to use distinct random values in different prime fields to reduce the risk of accidental cancellations.
The corresponding file `shift` may then read

```
443333693649116027 530304748736793718 4727957290386205426
2611273300125671761 4082324484597984495 3075032569815551545
...
```

The probes have then to be placed into the directory `probes`.
For each prime field, one requires a separate file `$PRIME.gz`, where `$PRIME` is the iterator number of the prime fields, i.e. 0, 1, 2, ...
The file `0.gz` for the first prime field may then look like

```
1 1 | 6303994061519960968 | 0 1 6747327755169076995
1 1 | 2956691294852519450 | 0 1 3400024988501635477
1 1 | 9142216511497419146 | 0 1 362178168291759390
2 1 | 8389701825437495663 | 0 1 8833035519086611690
1 2 | 1460015629339790903 | 0 1 1903349322988906930
...
```

where the three numbers in the first fold are the `zi_orders` (powers) `j_2`, `j_3`, `j_4` for the anchor points (`j_1` is irrelevant since the corresponding anchor point is 1).
The number in the second fold is `t_k`, while the numbers in the third fold are the results of the three black box functions in this example evaluated at the chosen values.

The executable can then be run with

```
firefly -v $NUMBER_OF_VARIABLES -p $THREADS
```

If the provided probes are not sufficient to fully interpolate and reconstruct the black-box functions, it will abort, store the intermediate results to the folder `ff_save`, and write the next probes it requires to the file `requested_probes.gz`, which for example reads

```
2 1 | 1931531046357483833 | 2374864740006599860 6216802671037139389 2374659357042260498
...
```

Again, the numbers in the first fold are the `zi_orders` (powers) of the anchor points.
The number of the second fold is a suggestion on which `t_k` to choose, while the numbers in the third fold are the values for the variables computed through the formula above assuming the suggested `t_k`.
The user is free to choose any other `t_k` if desired.
However, any `t_k` chosen twice for the same powers in the same prime field is automatically discarded by `FireFly`.

The user then has to evaluate the black box at those values and provide them in the `probes` directory as described above, e.g. by replacing the previous file.
By calling `firefly` again, the saved state is loaded from the `ff_save` directory and the interpolation resumed with the new probes.

In the first prime field, the functional forms of the black-box functions is not known.
Hence, `FireFly` only requests a relatively small number of probes after each step.
Once the functional forms are known in the second prime field, all probes required in the prime field are requested at the beginning of the prime field (after a first probe is fed).

A simple example can be found in the directory `examples/firefly_executable`.
