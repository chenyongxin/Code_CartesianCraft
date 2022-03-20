#[=======================================================================[.rst:
FindHYPRE
---------

Find HYPRE library.

This module will set the following variables in your project:

``HYPRE_FOUND``
  HYPRE was found on the system
``HYPRE_INCLUDE_DIR``
  HYPRE include directory
``HYPRE_LIBRARY``
  HYPRE library path

#]=======================================================================]

include(FindPackageHandleStandardArgs)

## Additional search paths
# Initialize variables
set(HYPRE_LIBRARY_PATH)
set(HYPRE_INCLUDE_PATH)

# Make sure the environment variable is not empty
# i.e. 4 arguments for string replacement.
if(NOT "$ENV{LIBRARY_PATH}" STREQUAL "")
	string(REPLACE ":" ";" HYPRE_LIBRARY_PATH $ENV{LIBRARY_PATH})
endif()
if(NOT "$ENV{CPATH}" STREQUAL "")
	string(REPLACE ":" ";" HYPRE_INCLUDE_PATH $ENV{CPATH})
endif()

find_path(HYPRE_INCLUDE_DIR
	NAMES HYPRE.h Hypre.h hypre.h
	HINTS
		/usr
		/usr/local
		$ENV{HYPRE_HOME}
		$ENV{HYPRE_DIR}
		${HYPRE_INCLUDE_PATH}
	PATH_SUFFIXES 
		include Include
)

find_library(HYPRE_LIBRARY
	NAMES 
		HYPRE Hypre hypre
	HINTS 
		/usr
		/usr/local
		$ENV{HYPRE}
		$ENV{HYPRE_HOME}
		$ENV{HYPRE_DIR}
		${HYPRE_LIBRARY_PATH}
	PATH_SUFFIXES 
		lib lib64 Lib
)

find_package_handle_standard_args(HYPRE REQUIRED_VARS HYPRE_LIBRARY HYPRE_INCLUDE_DIR)