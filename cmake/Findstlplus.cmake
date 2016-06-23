# - Tries to find LocCom library
# once found it defines:
# 
#  LIBSTLPLUS_FOUND - System has LibXml2
#  LIBSTLPLUS_INCLUDE_DIRS - The LibXml2 include directories
#  LIBSTLPLUS_LIBRARIES - The libraries needed to use LibXml2

find_path(STLPLUS_INCLUDE_DIR stlplus3/portability/portability.hpp PATH_SUFFIXES include)
find_library(STLPLUS_LIBRARY NAMES stlplus3_portabilty  PATH_SUFFIXES lib)

set(STLPLUS_INCLUDE_DIRS ${STLPLUS_INCLUDE_DIR})
set(STLPLUS_LIBRARIES ${STLPLUS_LIBRARY} ${CMAKE_DL_LIBS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibSTLPLUS  DEFAULT_MSG
	STLPLUS_INCLUDE_DIR STLPLUS_LIBRARY)

mark_as_advanced(STLPLUS_INCLUDE_DIR STLPLUS_LIBRARY )
