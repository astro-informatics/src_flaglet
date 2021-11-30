get_filename_component(Flaglet_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
message(STATUS "Linking to flaglet package in ${Flaglet_CMAKE_DIR}")
if(NOT TARGET flaglet AND EXISTS "${Flaglet_CMAKE_DIR}/FlagletTargets.cmake")
  include("${Flaglet_CMAKE_DIR}/FlagletTargets.cmake")
endif()

set(Flaglet_INCLUDE_DIRS "@ALL_INCLUDE_DIRS@")
set(Flaglet_LIBRARIES flaglet)
