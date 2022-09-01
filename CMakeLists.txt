cmake_minimum_required(VERSION 3.24)

project(numerics)
set(CMAKE_CXX_STANDARD 17)
file(GLOB SOURCES src/*)
add_library(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC include)

enable_testing()
add_executable(root_test tests/root_test.cpp)
target_link_libraries(root_test PRIVATE ${PROJECT_NAME})
add_test(NAME root_test COMMAND root_test)

add_executable(integral_test tests/integral_test.cpp)
target_link_libraries(integral_test PRIVATE ${PROJECT_NAME})
add_test(NAME integral_test COMMAND integral_test)