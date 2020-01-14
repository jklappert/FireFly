cmake_minimum_required(VERSION 3.1)
project(FireFly VERSION 1.3.4 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# decleare some colours
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

find_package(GMP REQUIRED)
find_package(ZLIB REQUIRED)

# set the build type to release. If one needs debugging symbols, compile
# with "-DCMAKE_BUILD_TYPE=Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable all warnings in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -pedantic -g")
endif()

# enable optimization
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -g")
endif()

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(INCLUDE_PATH ${SOURCE_PATH}/include)

# include headers
include_directories(${INCLUDE_PATH})
include_directories(${GMP_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS})
set(PKG_INCLUDE -I${GMP_INCLUDE_DIR} -I${ZLIB_INCLUDE_DIRS})

# Check for FLINT
# if you want to use FLINT instead of our default implementation, call cmake with
# cmake -DWITH_FLINT=true ..
# set(WITH_FLINT false)
if("${CUSTOM}" STREQUAL true)
  add_definitions(-DCUSTOM)
  set(PREPROCESSOR_FLAGS "-DCUSTOM")
elseif("${WITH_FLINT}" STREQUAL true)
  find_package(FLINT REQUIRED)
  include_directories(${FLINT_INCLUDE_DIR})
  set(PKG_INCLUDE ${PKG_INCLUDE} -I${FLINT_INCLUDE_DIR})
  add_definitions(-DFLINT)
  set(PREPROCESSOR_FLAGS "-DFLINT")
else()
  message("-- ${Red}Using default implementation of modular arithmetic. This can be much slower than FLINT!${ColourReset}")
  add_definitions(-DDEFAULT)
  set(PREPROCESSOR_FLAGS "-DDEFAULT")
endif()

# Check for MPI
if("${WITH_MPI}" STREQUAL true)
  if("${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
    find_package(MPI REQUIRED)
  else()
    message("-- Using user defined MPI implementation library: ${MPI_CXX_LIBRARIES}")
    message("-- Using user defined MPI implementation headers: ${MPI_CXX_INCLUDE_PATH}")
  endif()
  add_definitions(-DWITH_MPI)
  set(PREPROCESSOR_FLAGS "${PREPROCESSOR_FLAGS} -DWITH_MPI")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(PKG_INCLUDE ${PKG_INCLUDE} -I${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()

# Firefly library files
set(LIBSOURCES
  ${SOURCE_PATH}/BaseReconst.cpp
  ${SOURCE_PATH}/FFInt.cpp
  ${SOURCE_PATH}/FFIntVec.cpp
  ${SOURCE_PATH}/FFThieleInterpolator.cpp
  ${SOURCE_PATH}/gzstream.cpp
  ${SOURCE_PATH}/HornerGenerator.cpp
  ${SOURCE_PATH}/Monomial.cpp
  ${SOURCE_PATH}/ParserUtils.cpp
  ${SOURCE_PATH}/Poly.cpp
  ${SOURCE_PATH}/Polynomial.cpp
  ${SOURCE_PATH}/PolynomialFF.cpp
  ${SOURCE_PATH}/PolyReconst.cpp
  ${SOURCE_PATH}/RationalFunction.cpp
  ${SOURCE_PATH}/RationalFunctionFF.cpp
  ${SOURCE_PATH}/RationalNumber.cpp
  ${SOURCE_PATH}/RatReconst.cpp
  ${SOURCE_PATH}/ReconstHelper.cpp
  ${SOURCE_PATH}/ShuntingYardParser.cpp
  ${SOURCE_PATH}/utils.cpp)

set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/firefly")
set(INCLUDE_INSTALL_DIR_2 "${CMAKE_INSTALL_PREFIX}/include")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")

configure_file (
  "${INCLUDE_PATH}/version.hpp.in"
  "${INCLUDE_PATH}/version.hpp")

add_library(FireFly_static STATIC ${LIBSOURCES})
set_target_properties(FireFly_static PROPERTIES OUTPUT_NAME firefly)
add_library(FireFly_shared SHARED ${LIBSOURCES})
set_target_properties(FireFly_shared PROPERTIES OUTPUT_NAME firefly)

set(LIBS FireFly_static ${GMP_LIBRARIES} ${ZLIB_LIBRARIES})
if("${WITH_FLINT}" STREQUAL true)
  set(LIBS ${LIBS} ${FLINT_LIBRARY})
endif()

if("${WITH_MPI}" STREQUAL true)
  #string(STRIP ${MPI_CXX_LINK_FLAGS} MPI_CXX_LINK_FLAGS)
  #set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
  set(EXAMPLE_EXE ${PROJECT_SOURCE_DIR}/example_mpi.cpp)
else()
  set(EXAMPLE_EXE ${PROJECT_SOURCE_DIR}/example.cpp)
endif()

if("${CUSTOM}" STREQUAL "true")
  message("-- ${Red}Using custom modular arithmetic. The example is thus not compiled.${ColourReset}")
else()
  add_executable(example ${EXAMPLE_EXE})
  target_link_libraries(example ${LIBS})
endif()

# install libraries
install(
  TARGETS FireFly_static FireFly_shared
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

set(PKG_LIBS ${LIBS})
list(REMOVE_AT PKG_LIBS 0)
list(REMOVE_DUPLICATES PKG_LIBS)
string(REPLACE ";" " " PKG_LIBS "${PKG_LIBS}")

list(REMOVE_DUPLICATES PKG_INCLUDE)
string(REPLACE ";" " " PKG_INCLUDE "${PKG_INCLUDE}")

configure_file("cmake/pc.in" "firefly.pc" @ONLY)

# install pkgconfig file
install(FILES
  "${CMAKE_BINARY_DIR}/firefly.pc"

DESTINATION "${LIB_INSTALL_DIR}/pkgconfig")

# install public headers
install(FILES
  "${INCLUDE_PATH}/BaseReconst.hpp"
  "${INCLUDE_PATH}/BlackBoxBase.hpp"
  "${INCLUDE_PATH}/DenseSolver.hpp"
  "${INCLUDE_PATH}/FFInt.hpp"
  "${INCLUDE_PATH}/FFIntVec.hpp"
  "${INCLUDE_PATH}/FFThieleInterpolator.hpp"
  "${INCLUDE_PATH}/gzstream.hpp"
  "${INCLUDE_PATH}/HornerGenerator.hpp"
  "${INCLUDE_PATH}/Logger.hpp"
  "${INCLUDE_PATH}/Monomial.hpp"
  "${INCLUDE_PATH}/ParserUtils.hpp"
  "${INCLUDE_PATH}/Poly.hpp"
  "${INCLUDE_PATH}/Polynomial.hpp"
  "${INCLUDE_PATH}/PolynomialFF.hpp"
  "${INCLUDE_PATH}/PolyReconst.hpp"
  "${INCLUDE_PATH}/RationalFunction.hpp"
  "${INCLUDE_PATH}/RationalFunctionFF.hpp"
  "${INCLUDE_PATH}/RationalNumber.hpp"
  "${INCLUDE_PATH}/RatReconst.hpp"
  "${INCLUDE_PATH}/ReconstHelper.hpp"
  "${INCLUDE_PATH}/Reconstructor.hpp"
  "${INCLUDE_PATH}/ShuntingYardParser.hpp"
  "${INCLUDE_PATH}/ThreadPool.hpp"
  "${INCLUDE_PATH}/tinydir.h"
  "${INCLUDE_PATH}/UintHasher.hpp"
  "${INCLUDE_PATH}/utils.hpp"
  "${INCLUDE_PATH}/version.hpp"
  DESTINATION "${INCLUDE_INSTALL_DIR}")

if("${WITH_MPI}" STREQUAL "true")
  install(FILES "${INCLUDE_PATH}/MPIWorker.hpp" DESTINATION "${INCLUDE_INSTALL_DIR}")
endif()

# adding documentation
add_subdirectory(doc)

# test
#if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
#  enable_testing()
#  add_subdirectory(test)
#endif()
