# Exports flaglet so other packages can access it
export(TARGETS flaglet FILE "${PROJECT_BINARY_DIR}/FlagletTargets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE flaglet)
endif()

# First in binary dir
set(ALL_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
configure_File(cmake/FlagletConfig.in.cmake
  "${PROJECT_BINARY_DIR}/FlagletConfig.cmake" @ONLY
)
configure_File(cmake/FlagletConfigVersion.in.cmake
  "${PROJECT_BINARY_DIR}/FlagletConfigVersion.cmake" @ONLY
)

# Then for installation tree
file(RELATIVE_PATH REL_INCLUDE_DIR
    "${CMAKE_INSTALL_PREFIX}/share/cmake/flaglet"
    "${CMAKE_INSTALL_PREFIX}/include/flaglet"
)
set(ALL_INCLUDE_DIRS "\${Flaglet_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/FlagletConfig.in.cmake
  "${PROJECT_BINARY_DIR}/CMakeFiles/FlagletConfig.cmake" @ONLY
)

# Finally install all files
install(FILES
  "${PROJECT_BINARY_DIR}/CMakeFiles/FlagletConfig.cmake"
  "${PROJECT_BINARY_DIR}/FlagletConfigVersion.cmake"
    DESTINATION share/cmake/flaglet
    COMPONENT dev
)

install(EXPORT FlagletTargets DESTINATION share/cmake/flaglet COMPONENT dev)
