add_executable(tests main.cpp tests.cpp)
target_link_libraries(tests la_vector)
add_test(NAME test_all COMMAND tests)