# Try to find the standalone Ariadne Numeric library.
#
# Once done this will define:
#
#  AriadneNumeric_FOUND - system has Ariadne Numeric
#  ARIADNE_NUMERIC_INCLUDE_DIRS - the Ariadne Numeric include directories
#  ARIADNE_NUMERIC_LIBRARIES - libraries required to use Ariadne Numeric
#
# and, when found:
#
#  AriadneNumeric::ariadne-numeric - imported library target

# This file is part of Ariadne.

# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.

find_library(ARIADNE_NUMERIC_LIBRARY NAMES ariadne-numeric)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(GMP QUIET gmp)
  pkg_check_modules(MPFR QUIET mpfr)
endif()

find_path(GMP_INCLUDE_DIR gmp.h HINTS ${GMP_INCLUDE_DIRS})
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${MPFR_INCLUDE_DIRS})
find_library(GMP_LIBRARY NAMES gmp libgmp HINTS ${GMP_LIBRARY_DIRS})
find_library(MPFR_LIBRARY NAMES mpfr libmpfr HINTS ${MPFR_LIBRARY_DIRS})

find_path(ARIADNE_NUMERIC_INCLUDE_DIR ariadne-numeric.hpp PATH_SUFFIXES ariadne-numeric)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneNumeric DEFAULT_MSG
  ARIADNE_NUMERIC_LIBRARY
  ARIADNE_NUMERIC_INCLUDE_DIR
  GMP_INCLUDE_DIR
  MPFR_INCLUDE_DIR
  GMP_LIBRARY
  MPFR_LIBRARY
)

set(ARIADNE_NUMERIC_FOUND ${AriadneNumeric_FOUND})

if(AriadneNumeric_FOUND)
  get_filename_component(ARIADNE_NUMERIC_INCLUDE_PARENT_DIR ${ARIADNE_NUMERIC_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_NUMERIC_INCLUDE_DIRS
    ${ARIADNE_NUMERIC_INCLUDE_PARENT_DIR}
    ${ARIADNE_NUMERIC_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
  )
  set(ARIADNE_NUMERIC_LIBRARIES
    ${ARIADNE_NUMERIC_LIBRARY}
    ${MPFR_LIBRARY}
    ${GMP_LIBRARY}
  )

  if(NOT TARGET AriadneNumeric::ariadne-numeric)
    add_library(AriadneNumeric::ariadne-numeric UNKNOWN IMPORTED)
    set_target_properties(AriadneNumeric::ariadne-numeric PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_NUMERIC_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_NUMERIC_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${MPFR_LIBRARY};${GMP_LIBRARY}"
    )
  endif()
endif()

mark_as_advanced(
  ARIADNE_NUMERIC_INCLUDE_DIR
  ARIADNE_NUMERIC_LIBRARY
  GMP_INCLUDE_DIR
  GMP_LIBRARY
  MPFR_INCLUDE_DIR
  MPFR_LIBRARY
)
