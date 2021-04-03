cmake_minimum_required(VERSION 3.5)
project(Global_Electric_Circuit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O2")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -O2")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

enable_language(Fortran)

add_subdirectory(irbem)

add_library(cnpy SHARED "cnpy/cnpy.cpp")

add_executable(Global_Electric_Circuit
        main.cpp
        integral.cpp
        integral.h
        std_atm.h
        sigma.h
        geomagnetic.cpp
        geomagnetic.h)

target_link_libraries(Global_Electric_Circuit
        PUBLIC ${Boost_LIBRARIES}
        irbem
        cnpy
        )
target_include_directories(Global_Electric_Circuit SYSTEM PUBLIC ${Boost_INCLUDE_DIR})
