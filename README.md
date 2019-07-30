# FireFly

[![Build Status](https://travis-ci.org/jklappert/FireFly.svg?branch=master)](https://travis-ci.org/jklappert/FireFly) [![](https://img.shields.io/github/tag/jklappert/firefly)](https://gitlab.com/firefly-library/firefly/-/tags/1.2.0)

FireFly is a reconstruction library for rational functions and polynomials written in C++.

Please refer to this paper when using FireFly:
* J. Klappert and F. Lange, *Reconstructing Rational Functions with FireFly*, [[1904.00009](https://arxiv.org/abs/1904.00009)]

## Requirements
FireFly requires:
* C++ compiler supporting C++11
* [CMake](https://cmake.org/) >= 3.1
* [FLINT](http://www.flintlib.org/) >= 2.5 (optional)
* [GMP](https://gmplib.org/) >= 6.1

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

By default the code is compiled with optimizations.

If FLINT is used for modular arithmetic and it cannot be found in the default system directories, one has to add the additional flags:

```
-DFLINT_INCLUDE_DIR=$FLINT_INC_PATH -DFLINT_LIBRARY=$FLINT_LIB_PATH
```

where `FLINT_LIB_PATH` is the absolute path pointing to the shared library of FLINT.


## Reconstructing functions
To reconstruct functions with FireFly it offers an interface which directly makes use of a thread pool for the parallel reconstruction of various functions over the same prime field. Additionally, black-box probes are calculated in parallel.

The black box is implemented as functor. The user has to define the black box as a derived class of `BlackBoxBase` and provide a constructor, the evaluation of the black box, and a function which allows the user to change member variables when the prime field changes:

```cpp
class BlackBoxUser : public BlackBoxBase {
public:
  BlackBoxUser(...) {
    ...
  }
  virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values) {
    ...
  }
  virtual void prime_changed() {
    ...
  }
}
```

The actual reconstruction is done by the `Reconstructor` object which is initialized with the number of variables `n_var`, the number of threads `n_thr` that should be used, and the `BlackBoxUser` object `bb`:

```cpp
BlackBoxUser bb(...);
Reconstructor rec(n_var, n_thr, bb);
```

The reconstruction is then started by

```cpp
rec.reconstruct();
```

The reconstruction will run from this point until it is finished. The results can be obtained with

```cpp
rec.get_result();
```

Additional options can be set and we refer to the `example.cpp` file and the code documentation.


## Directly parse collections of rational functions
When the conversion and compilation steps of the `convert_to_ff.sh` scripts are too time consuming, FireFly also provides a parser class for rational functions. The functions will be stored in reverse polish notation to be evaluated for a given parameter point. Parsing collections of rational functions can be done with the `ShuntingYardParser` class. It has to be constructed with a path to a file which contains the rational functions and a vector which sets the occurring variables:

```cpp
ShuntingYardParser parser(path, vars);
```

Here, `path` is a string containing the path to the file in which the needed functions are stored and `vars` is a vector of strings which represent the occurring variables. The collection of functions have to be separated by new lines to be identified correctly, e.g.,

```
(x+y*3)/(z^2+x)
(12132132323213213212/33*x + 12 * 3 - x^100*y^2)/(3*y^5)
```

The corresponding vector `vars` would thus be

```
std::vector<std::string> vars = {"x","y","z"};
```

The functions have to be parsed only once and can be evaluated afterwards calling

```cpp
parser.evaluate(values);
//parser.evaluate_pre(values); // Evaluates the black-box functions with precomputed values (faster than evaluate()). Requires parser.precompute_tokens() after the field has changed.
```

where `values` is a vector which contains the parameter point at which the functions should be evaluated. The function `evaluate` returns a vector of `FFInt` objects which is filled by the values of the evaluated functions in the same order as the functions are defined in the input file. Thus, it can be directly used in the `BlackBox` functor of FireFly. An example file is given in `s_y_test.m`. Note that only the operators

```
+, -, *, /, ^
```

are supported.

For convenience, FireFly also provides a script which converts a list of rational functions (stored as an expression list of Mathematica) to FireFly's parsable format. It is located in the `mma_2_ff` directory and can be executed with

```
./convert_to_sy.sh $FILE
```

where `$FILE` contains the list of rational functions.


## Converting Mathematica expressions to C++ code
Sometimes the black box is not provided by a code but some Mathematica expressions. For this purpose FireFly provides a script to convert Mathematica functions to compilable C++ code. This can be useful by performing algebraic computations on large functions (see also the parser of FireFly). The functions have to be provided as a file in which a list of functions (expression or string) is stored, e.g.,

```
{x+y,2*x+z,...}
```

Additionally, a file is needed in which the occurring variables are stored as a Mathematica list, e.g.,

```
{x,y,z,...}
```

Note that both lists are allowed to contain expressions and/or strings. The scripts are located in the `mma_2_ff` directory. One can than generate C++ code by calling

```
cd mma_2_ff
chmod u+x convert_to_ff.sh
./convert_to_ff.sh $PATH_TO_FUNCTION_FILE $PATH_TO_VARIABLES_FILE <number_of_threads>
```

Here, `$PATH_TO_FUNCTION_FILE` defines the path to the file in which the list of functions is stored, `$PATH_TO_VARIABLES_FILE` defines the path to the file in which the list of variables is stored, and `<number_of_threads>` defines the number of threads you want to use for the reconstruction process. The latter is given as an integer. During the conversion process a directory `ff_conv` is created in which the C++ files are written. It contains an executable `exec.cpp` which has to be modified to your needs (filling the black box, doing something with the result,...) and all functions which are splitted to numerator and denominator and to subfunctions if they exceed a length of 500 terms. The functions are numbered according to their ordering in the list and can be evaluated in C++ by calling, e.g.,

```cpp
std::vector<FFInt> values = {1, 5, 7};
FFInt res = fun1(values) + fun2(values) + ...;
```

When using the `black_box` functor in `exec.cpp`, the vector `values` will be filled by FireFly. After a specification what the black box should be, the generated makefile will compile your program by calling

```
make
```

A build directory will be created (`/build`) in which the executable and the object files can be found. To execute the reconstruction one just has to call

```
./build/exec
```

which will be performed using the number of threads defined in the conversion process.


## Code Documentation
Doxygen can be used to generate code documentation. To generate the documentation, run

```
make doc
```

The generated documentation can be found in `doc/html/index.html`.
