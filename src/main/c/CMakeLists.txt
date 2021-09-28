add_library(flaglet STATIC flaglet_axisym.c flaglet_tiling.c flaglet_transform.c)
target_link_libraries(flaglet PUBLIC s2let::s2let GSL::gsl GSL::gslcblas ${FLAG_TARGET})
target_include_directories(
  flaglet
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<INSTALL_INTERFACE:include> ${Flag_INCLUDE_DIR})
set_target_properties(flaglet PROPERTIES C_STANDARD 99 ARCHIVE_OUTPUT_DIRECTORY
                                                   ${PROJECT_BINARY_DIR}/lib)
configure_file(${PROJECT_SOURCE_DIR}/include/flaglet_version.in.h
               ${PROJECT_BINARY_DIR}/include/flaglet_version.h)

if(fPIC)
  set_target_properties(flaglet PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

if(NOT SKBUILD)
  install(
    TARGETS flaglet
    EXPORT FlagletTargets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER)

  install(
    FILES ${PROJECT_SOURCE_DIR}/include/flaglet.h
          ${PROJECT_SOURCE_DIR}/include/flaglet_axisym.h
          ${PROJECT_SOURCE_DIR}/include/flaglet_tiling.h
          ${PROJECT_SOURCE_DIR}/include/flaglet_transform.h
          ${PROJECT_SOURCE_DIR}/include/flaglet_types.h
    DESTINATION include/flaglet)
endif()
