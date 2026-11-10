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

find_path(ARIADNE_NUMERIC_INCLUDE_DIR ariadne-numeric.hpp PATH_SUFFIXES ariadne-numeric)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneNumeric DEFAULT_MSG
  ARIADNE_NUMERIC_LIBRARY
  ARIADNE_NUMERIC_INCLUDE_DIR
)

set(ARIADNE_NUMERIC_FOUND ${AriadneNumeric_FOUND})

if(AriadneNumeric_FOUND)
  get_filename_component(ARIADNE_NUMERIC_INCLUDE_PARENT_DIR ${ARIADNE_NUMERIC_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_NUMERIC_INCLUDE_DIRS
    ${ARIADNE_NUMERIC_INCLUDE_PARENT_DIR}
    ${ARIADNE_NUMERIC_INCLUDE_DIR}
  )
  set(ARIADNE_NUMERIC_LIBRARIES
    ${ARIADNE_NUMERIC_LIBRARY}
  )

  if(NOT TARGET AriadneNumeric::ariadne-numeric)
    add_library(AriadneNumeric::ariadne-numeric UNKNOWN IMPORTED)
    set_target_properties(AriadneNumeric::ariadne-numeric PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_NUMERIC_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_NUMERIC_INCLUDE_DIRS}"
    )
  endif()
endif()

mark_as_advanced(
  ARIADNE_NUMERIC_INCLUDE_DIR
  ARIADNE_NUMERIC_LIBRARY
)
