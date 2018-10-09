cmake_minimum_required(VERSION 3.1)
project(Firefly VERSION 0.0.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

find_package(GMP REQUIRED)

# set your favorite compile here
set (CMAKE_CXX_COMPILER clang++)

# to use Valgrind uncomment the following line
#add_definitions(-DVALGRIND -g)


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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")
  endif()
endif()

# enable optimization
if(CMAKE_BUILD_TYPE STREQUAL Release)
  if(CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")
  endif()
endif()

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(INCLUDE_PATH ${SOURCE_PATH}/include)

# include headers
include_directories(${INCLUDE_PATH})

# Firefly library files
file(GLOB LIBSOURCES
  ${SOURCE_PATH}/*.cpp)

add_library(Firefly STATIC ${LIBSOURCES})

add_executable(example ${PROJECT_SOURCE_DIR}/example.cpp)
target_link_libraries(example Firefly gmp)

# needed to install libray
#install(TARGETS Firefly_static RUNTIME DESTINATION bin)
