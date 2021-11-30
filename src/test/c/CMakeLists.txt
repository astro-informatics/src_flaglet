foreach(testname flaglet)
  add_executable(test_${testname} test_${testname}.c)
  target_link_libraries(test_${testname} PRIVATE flaglet cmocka)
  set_target_properties(
    test_${testname} PROPERTIES C_STANDARD 11 RUNTIME_OUTPUT_DIRECTORY
                                              ${PROJECT_BINARY_DIR}/bin)
  add_test(NAME test_${testname} COMMAND test_${testname})
endforeach()
