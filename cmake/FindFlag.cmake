find_path(Flag_INCLUDE_DIR NAMES flag/flag.h)
find_library(Flag_LIBRARY NAMES flag)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Flag  DEFAULT_MSG Flag_INCLUDE_DIR Flag_LIBRARY)
MESSAGE(STATUS "Using flag include ${Flag_INCLUDE_DIR}")
MESSAGE(STATUS "Using flag lib ${Flag_LIBRARY}")
if(NOT FLAG_FOUND)
  set(Flag_INCLUDE_DIR "")
  set(Flag_LIBRARY "")
endif()
