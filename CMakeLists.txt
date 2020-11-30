cmake_minimum_required(VERSION 3.5)

project(Global_Electric_Circuit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Global_Electric_Circuit main.cpp integral.cpp integral.h StdAtm.h sigma.h)

target_link_libraries(Global_Electric_Circuit PUBLIC ${Boost_LIBRARIES})
target_include_directories(Global_Electric_Circuit SYSTEM PUBLIC ${Boost_INCLUDE_DIR})
