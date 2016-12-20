cmake_minimum_required(VERSION 3.6)
project(EigenTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(SOURCE_FILES main.cpp MDsystem.h MDintegrator.h VelocityVerletIntegrator.cpp VelocityVerletIntegrator.h MDtracker.cpp MDtracker.h Argon.cpp Argon.h Histogram.cpp Histogram.h)
add_executable(EigenTest ${SOURCE_FILES})
add_executable(main main.cpp MDsystem.h MDintegrator.h VelocityVerletIntegrator.cpp VelocityVerletIntegrator.h MDtracker.cpp MDtracker.h Argon.cpp Argon.h Histogram.cpp Histogram.h)
