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
* [Parse collections of rational functions](#parse-collections-of-rational-functions)
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
  The functions as well as the source code can be found in the directory `benchmarks`.
  Both executables can simply be run without command line arguments.
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


## Parse collections of rational functions
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

We refer to the [Wiki page](https://gitlab.com/firefly-library/firefly/-/wikis/FireFly-Executable-(in-Development)) for detailed instructions.
