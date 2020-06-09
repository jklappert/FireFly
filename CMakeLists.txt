cmake_minimum_required(VERSION 3.1)
project(FireFly VERSION 2.0.1 LANGUAGES CXX)

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

find_package(GMP 6.1.2 REQUIRED)
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -pthread -g")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -g")
endif()

set(CMAKE_SHARED_LINKER_FLAGS "-pthread")

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
add_definitions("-DFIREFLY_ROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
set(INCLUDE_PATH ${SOURCE_PATH}/include)
set(INCLUDE_PATH_SRC ${SOURCE_PATH}/include/firefly)

# include headers
include_directories(${INCLUDE_PATH})
include_directories(${GMP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
set(PKG_INCLUDE -I${GMP_INCLUDE_DIRS} -I${ZLIB_INCLUDE_DIRS})

set(FLINT_FLAG "0")
set(DEFAULT_FLAG "0")
set(CUSTOM_FLAG "0")
set(MPI_FLAG "0")

# Check for FLINT
# if you want to use FLINT instead of our default implementation, call cmake with
# cmake -DWITH_FLINT=true ..
# set(WITH_FLINT false)
if("${CUSTOM}" STREQUAL true)
#  add_definitions(-DCUSTOM)
  set(PREPROCESSOR_FLAGS "-DCUSTOM")
  set(CUSTOM_FLAG "1")
elseif("${WITH_FLINT}" STREQUAL true)
  find_package(FLINT REQUIRED)
  include_directories(${FLINT_INCLUDE_DIR})
  set(PKG_INCLUDE ${PKG_INCLUDE} -I${FLINT_INCLUDE_DIR})
#  add_definitions(-DFLINT)
  set(PREPROCESSOR_FLAGS "-DFLINT")
  set(FLINT_FLAG "1")
else()
  message("-- ${Red}Using default implementation of modular arithmetic. This can be much slower than FLINT!${ColourReset}")
#  add_definitions(-DDEFAULT)
  set(PREPROCESSOR_FLAGS "-DDEFAULT")
  set(DEFAULT_FLAG "1")
endif()

# Check for MPI
if("${WITH_MPI}" STREQUAL true)
  set(MPI_FLAG "1")
  if("${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
    find_package(MPI REQUIRED)
  else()
    message("-- Using user-defined MPI implementation library: ${MPI_CXX_LIBRARIES}")
    message("-- Using user-defined MPI implementation headers: ${MPI_CXX_INCLUDE_PATH}")
  endif()
#  add_definitions(-DWITH_MPI)
  set(PREPROCESSOR_FLAGS "${PREPROCESSOR_FLAGS} -DWITH_MPI")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(PKG_INCLUDE ${PKG_INCLUDE} -I${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()

# Firefly library files
set(LIBSOURCES
  ${SOURCE_PATH}/AmplitudeParser.cpp
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
set(BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INCLUDE_INSTALL_DIR_2 "${CMAKE_INSTALL_PREFIX}/include")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")

configure_file (
  "${INCLUDE_PATH_SRC}/version.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/version.hpp")

configure_file (
  "${INCLUDE_PATH_SRC}/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/config.hpp")

add_library(FireFly_static STATIC ${LIBSOURCES})
set_target_properties(FireFly_static PROPERTIES OUTPUT_NAME firefly)
add_library(FireFly_shared SHARED ${LIBSOURCES})
set_target_properties(FireFly_shared PROPERTIES OUTPUT_NAME firefly)

set(LIBS ${GMP_LIBRARIES} ${ZLIB_LIBRARIES})
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
  set(INSERT_IBPS_EXE ${PROJECT_SOURCE_DIR}/ff_insert.cpp)
endif()

if("${CUSTOM}" STREQUAL "true")
  message("-- ${Red}Using custom modular arithmetic. The example is thus not compiled.${ColourReset}")
elseif("${WITH_MPI}" STREQUAL "true")
  add_executable(example ${EXAMPLE_EXE})
  if ("${WITH_JEMALLOC}" STREQUAL "true")
    execute_process(COMMAND "jemalloc-config" "--libdir" OUTPUT_VARIABLE JEMALLOC_LIB_DIR)
    execute_process(COMMAND "jemalloc-config" "--libs" OUTPUT_VARIABLE JEMALLOC_LIBS)
    string(REGEX REPLACE "\n$" "" JEMALLOC_LIB_DIR "${JEMALLOC_LIB_DIR}")
    string(REGEX REPLACE "\n$" "" JEMALLOC_LIBS "${JEMALLOC_LIBS}")
    message("-- Found jemalloc: ${JEMALLOC_LIB_DIR}")
    set(JEMALLOC_LIB "-L${JEMALLOC_LIB_DIR} -Wl,-rpath,${JEMALLOC_LIB_DIR} -ljemalloc" ${JEMALLOC_LIBS})
    set(LIBS ${JEMALLOC_LIB} ${LIBS})
  endif()
  target_link_libraries(example FireFly_static ${LIBS})
  target_link_libraries(FireFly_shared ${LIBS})
else()
  add_executable(example ${EXAMPLE_EXE})
  add_executable(ff_insert ${INSERT_IBPS_EXE})
  if ("${WITH_JEMALLOC}" STREQUAL "true")
    execute_process(COMMAND "jemalloc-config" "--libdir" OUTPUT_VARIABLE JEMALLOC_LIB_DIR)
    execute_process(COMMAND "jemalloc-config" "--libs" OUTPUT_VARIABLE JEMALLOC_LIBS)
    string(REGEX REPLACE "\n$" "" JEMALLOC_LIB_DIR "${JEMALLOC_LIB_DIR}")
    string(REGEX REPLACE "\n$" "" JEMALLOC_LIBS "${JEMALLOC_LIBS}")
    message("-- Found jemalloc: ${JEMALLOC_LIB_DIR}")
    set(JEMALLOC_LIB "-L${JEMALLOC_LIB_DIR} -Wl,-rpath,${JEMALLOC_LIB_DIR} -ljemalloc" ${JEMALLOC_LIBS})
    set(LIBS ${JEMALLOC_LIB} ${LIBS})
  endif()
  target_link_libraries(ff_insert FireFly_static ${LIBS})
  target_link_libraries(example FireFly_static ${LIBS})
  target_link_libraries(FireFly_shared ${LIBS})
endif()

if("${DISABLE_INSTALL}" STREQUAL "")
  # install libraries
  install(
    TARGETS FireFly_static FireFly_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)

  set(PKG_LIBS ${LIBS})
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
    "${INCLUDE_PATH_SRC}/AmplitudeParser.hpp"
    "${INCLUDE_PATH_SRC}/BaseReconst.hpp"
    "${INCLUDE_PATH_SRC}/BlackBoxBase.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/config.hpp"
    "${INCLUDE_PATH_SRC}/DenseSolver.hpp"
    "${INCLUDE_PATH_SRC}/FFInt.hpp"
    "${INCLUDE_PATH_SRC}/FFIntVec.hpp"
    "${INCLUDE_PATH_SRC}/FFThieleInterpolator.hpp"
    "${INCLUDE_PATH_SRC}/gzstream.hpp"
    "${INCLUDE_PATH_SRC}/HornerGenerator.hpp"
    "${INCLUDE_PATH_SRC}/Logger.hpp"
    "${INCLUDE_PATH_SRC}/Monomial.hpp"
    "${INCLUDE_PATH_SRC}/ParserUtils.hpp"
    "${INCLUDE_PATH_SRC}/Poly.hpp"
    "${INCLUDE_PATH_SRC}/Polynomial.hpp"
    "${INCLUDE_PATH_SRC}/PolynomialFF.hpp"
    "${INCLUDE_PATH_SRC}/PolyReconst.hpp"
    "${INCLUDE_PATH_SRC}/RationalFunction.hpp"
    "${INCLUDE_PATH_SRC}/RationalFunctionFF.hpp"
    "${INCLUDE_PATH_SRC}/RationalNumber.hpp"
    "${INCLUDE_PATH_SRC}/RatReconst.hpp"
    "${INCLUDE_PATH_SRC}/ReconstHelper.hpp"
    "${INCLUDE_PATH_SRC}/Reconstructor.hpp"
    "${INCLUDE_PATH_SRC}/ShuntingYardParser.hpp"
    "${INCLUDE_PATH_SRC}/ThreadPool.hpp"
    "${INCLUDE_PATH_SRC}/tinydir.h"
    "${INCLUDE_PATH_SRC}/UintHasher.hpp"
    "${INCLUDE_PATH_SRC}/utils.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/version.hpp"
    DESTINATION "${INCLUDE_INSTALL_DIR}")

  # install ibp insertion executable
  if(NOT "${CUSTOM}" STREQUAL "true" AND NOT "${WITH_MPI}" STREQUAL "true")
    install(TARGETS ff_insert DESTINATION "${BIN_INSTALL_DIR}")
  endif()

  if("${WITH_MPI}" STREQUAL "true")
    install(FILES "${INCLUDE_PATH_SRC}/MPIWorker.hpp" DESTINATION "${INCLUDE_INSTALL_DIR}")
  endif()

  # adding documentation
  add_subdirectory(doc)

  # test
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
    enable_testing()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
    add_subdirectory(test)
  endif()

  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To run the example type: ./example")
  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To run the inserter type: ./ff_insert")
  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
  add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To compile and run tests type: make check")
endif()
