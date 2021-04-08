cmake_minimum_required(VERSION 3.12)

project(
  flaglet
  VERSION "0.0.1"
  DESCRIPTION "Fast wavelets on the ball"
  HOMEPAGE_URL "https://github.com/astro-informatics/flaglet"
  LANGUAGES C)

option(tests "Enable testing" ON)
option(cfitsio "Links to cfitsio" OFF)
option(fPIC "Compile with fPIC" OFF)
option(conan_deps "Download dependencies using conan" ON)

find_package(GSL REQUIRED)
find_package(Flag REQUIRED)
set(FLAG_TARGET flag)

if(WIN32)
  set(fpic OFF)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include("sanitizers")
if(conan_deps OR CONAN_EDITABLE_MODE)
  include("conan_dependencies")
  if(NOT CONAN_EDITABLE_MODE)
    set(S2LET_TARGET CONAN_PKG::s2let)
  endif()
endif()
if(NOT conan_deps)
  find_package(S2let REQUIRED)
  set(S2LET_TARGET s2let)
endif()

add_subdirectory(src/main/c)

if(tests)
  enable_testing()
  include("fetch_cmocka")
  add_subdirectory(src/test/c)
endif()

if(NOT SKBUILD)
  include("exporting")
else()
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/main/pyflaglet)
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    set(DOXYGEN_FILE_PATTERNS *.h *.c)
    doxygen_add_docs(docs src/main/c src/test/c include)
  endif()
endif()
