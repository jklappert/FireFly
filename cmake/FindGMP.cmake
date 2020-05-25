# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
# Copyright (c) 2007, Francesco Biscani, <bluescarni@gmail.com>
# Copyright (c) 2019, Jonas Klappert and Fabian Lange

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ------------------------------------------------------------------------------------------

# Try to find the GMP libraries used with C++ code:
# GMP_FOUND - System has GMP lib
# GMP_INCLUDE_DIRS - The GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP
#
# Reads the following variables:
#
# - GMP_INCLUDE_DIR : include directory to search in
# - GMP_LIBRARY     : library directory to search in
#
# Sets the following variables:
#
# - GMP_FOUND       : System has GMP library
# - GMP_INCLUDE_DIRS : The GMP include directory
# - GMP_LIBRARIES   : Libraries needed to use GMP

# macro to extract the GMP version from a header file
macro(_FF_extract_GMP_version header_file)
  file(READ "${header_file}" _gmp_version_header)

  string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _gmp_version_major_match "${_gmp_version_header}")
  set(GMP_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)" _gmp_version_minor_match "${_gmp_version_header}")
  set(GMP_VERSION_MINOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _gmp_version_patch_match "${_gmp_version_header}")
  set(GMP_VERSION_RELEASE "${CMAKE_MATCH_1}")

  set(GMP_VERSION ${GMP_VERSION_MAJOR}.${GMP_VERSION_MINOR}.${GMP_VERSION_RELEASE})

  if(GMP_FIND_VERSION_EXACT AND NOT ${GMP_VERSION} VERSION_EQUAL ${GMP_FIND_VERSION})
    message(FATAL_ERROR "GMP version ${GMP_VERSION} found in ${GMP_H_INCLUDE_DIRS}, "
      "but exact version ${GMP_FIND_VERSION} is required.")
  elseif(${GMP_VERSION} VERSION_LESS ${GMP_FIND_VERSION})
    message(FATAL_ERROR "GMP version ${GMP_VERSION} found in ${GMP_H_INCLUDE_DIRS}, "
      "but at least version ${GMP_FIND_VERSION} is required.")
  endif()
endmacro()

# find version if asked for
if(GMP_FIND_VERSION AND GMP_INCLUDE_DIRS)
  # search for gmp.h
  find_path(GMP_H_INCLUDE_DIRS
    NAMES gmp.h
    PATHS ${GMP_INCLUDE_DIR}
    PATH_SUFFIXES include
  )

  if(GMP_H_INCLUDE_DIRS)
    _FF_extract_GMP_version("${GMP_H_INCLUDE_DIRS}/gmp.h")
  else()
    message(FATAL_ERROR "GMP version header gmp.h not found.")
  endif()

MESSAGE("-- Set GMP include path to: ${GMP_INCLUDE_DIRS}")
else()
  find_path(GMP_INCLUDE_DIRS
    NAMES gmpxx.h
    PATHS ${GMP_INCLUDE_DIR}
    PATH_SUFFIXES include
  )

  find_library(GMP_LIBRARIES
    NAMES gmp
    PATHS ${GMP_LIBRARY}
  )

  # find version if asked for
  if(GMP_FIND_VERSION AND GMP_INCLUDE_DIRS)
    # search for gmp.h
    find_path(GMP_H_INCLUDE_DIRS
      NAMES gmp.h
      PATHS ${GMP_INCLUDE_DIR}
      PATH_SUFFIXES include
    )
    if(GMP_H_INCLUDE_DIRS)
      _FF_extract_GMP_version("${GMP_H_INCLUDE_DIRS}/gmp.h")
    else()
      message(FATAL_ERROR "GMP version header gmp.h not found.")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP
  FOUND_VAR GMP_FOUND
  REQUIRED_VARS
    GMP_LIBRARIES
    GMP_INCLUDE_DIRS
)

mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARIES)
