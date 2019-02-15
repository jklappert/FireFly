cmake_minimum_required(VERSION 3.1)
project(FireFly VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

find_package(GMP REQUIRED)
find_package(FLINT REQUIRED)

# set your favorite compile here
set (CMAKE_CXX_COMPILER g++)

set(CMAKE_BUILD_TYPE Release)
# set the build type to release. If one needs debugging symbols, compile
# with "CMAKE_BUILD_TYPE Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable all warnings in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  if(CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -pedantic -g")
  endif()
endif()

# enable optimization
if(CMAKE_BUILD_TYPE STREQUAL Release)
  if(CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -g")
  endif()
endif()

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(INCLUDE_PATH ${SOURCE_PATH}/include)

# include headers
include_directories(${INCLUDE_PATH})

# Firefly library files
file(GLOB LIBSOURCES
  ${SOURCE_PATH}/BaseReconst.cpp
  ${SOURCE_PATH}/FFInt.cpp
  ${SOURCE_PATH}/Monomial.cpp
  ${SOURCE_PATH}/Polynomial.cpp
  ${SOURCE_PATH}/PolynomialFF.cpp
  ${SOURCE_PATH}/PolyReconst.cpp
  ${SOURCE_PATH}/RationalFunction.cpp
  ${SOURCE_PATH}/RationalNumber.cpp
  ${SOURCE_PATH}/RatReconst.cpp
  ${SOURCE_PATH}/ReconstHelper.cpp
  ${SOURCE_PATH}/Tests.cpp
  ${SOURCE_PATH}/utils.cpp)

file(GLOB LIBSOURCES_SHRD
  ${SOURCE_PATH}/BaseReconst.cpp
  ${SOURCE_PATH}/FFInt.cpp
  ${SOURCE_PATH}/Monomial.cpp
  ${SOURCE_PATH}/Polynomial.cpp
  ${SOURCE_PATH}/PolynomialFF.cpp
  ${SOURCE_PATH}/PolyReconst.cpp
  ${SOURCE_PATH}/RationalFunction.cpp
  ${SOURCE_PATH}/RationalNumber.cpp
  ${SOURCE_PATH}/RatReconst.cpp
  ${SOURCE_PATH}/ReconstHelper.cpp
  ${SOURCE_PATH}/utils.cpp)

set(EXEC_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/firefly")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")

add_library(FireFly_static STATIC ${LIBSOURCES})
set_target_properties(FireFly_static PROPERTIES OUTPUT_NAME firefly)
#add_library(FireFly_shared SHARED ${LIBSOURCES_SHRD})
#set_target_properties(FireFly_shared PROPERTIES OUTPUT_NAME firefly)

add_executable(example ${PROJECT_SOURCE_DIR}/example.cpp)
target_link_libraries(example FireFly_static gmp flint)

# install libraries
install(
  TARGETS FireFly_static# FireFly_shared
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

# install public headers
install(FILES
  "${INCLUDE_PATH}/BaseReconst.hpp"
  "${INCLUDE_PATH}/FFInt.hpp"
  "${INCLUDE_PATH}/Logger.hpp"
  "${INCLUDE_PATH}/Monomial.hpp"
  "${INCLUDE_PATH}/Polynomial.hpp"
  "${INCLUDE_PATH}/PolynomialFF.hpp"
  "${INCLUDE_PATH}/PolyReconst.hpp"
  "${INCLUDE_PATH}/RationalFunction.hpp"
  "${INCLUDE_PATH}/RationalNumber.hpp"
  "${INCLUDE_PATH}/RatReconst.hpp"
  "${INCLUDE_PATH}/ReconstHelper.hpp"
  "${INCLUDE_PATH}/UintHasher.hpp"
  "${INCLUDE_PATH}/utils.hpp"
  DESTINATION "${INCLUDE_INSTALL_DIR}")
