# FireFly

[![Build Status](https://travis-ci.org/jklappert/FireFly.svg?branch=master)](https://travis-ci.org/jklappert/FireFly)

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
FireFly uses CMake to generate files for build automation. To build FireFly one should first create a seperate `build` directory inside FireFly's top directory. Afterwards, `cmake` should be called:
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
To reconstruct functions with FireFly it offers an interface which directly makes use of a thread pool for the parallel reconstruction of various functions over the same prime field. Additionaly, black-box probes are calculated parallelized.

The reconstruction starts with initializing a `Reconstructor` object

```cpp
Reconstructor rec(n_var, n_thr);
```

with `n_var` being the number of variables and `n_thr` are the number of threads that should be used. Since the black-box function is called inside the `Reconstructor` class, it should be defined as

```cpp
void Reconstructor::black_box(vector<FFInt> result, const vector<FFInt>& values){
    ...
}
```

The `FFInt` object is used for all modular arithmetic operations. The STL vector `result` represents the black-box probes for a given tuple of values which are given in the STL vector `values`. It has to provide an immutable ordering and is filled by the user. After defining the black-box function, one could start the reconstruction by 

```cpp
rec.reconstruct();
```

The reconstruction will run from this point until it is finished. Additional options can be set and we refer to the `example.cpp` file and the code documentation.

## Code Documentation
Doxygen can be used to generate code documentation. To generate the documentation, run
```
make doc
```

The generated documentation can be found in `doc/html/index.html`.
