cmake_minimum_required (VERSION 3.5.1)
project (qcasim)

#system specific setting
if(     CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message("OS is Linux!!!")

elseif( CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message("OS is Mac OS X!!!")

elseif( CMAKE_SYSTEM_NAME STREQUAL "Windows")
    message("OS is Windows!!!")
endif()

#directory setting
include_directories($ENV{COMMON}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)

#build flags setting
set(CMAKE_CXX_FLAGS         "$ENV{CXXFLAGS} -Wall -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -O0 -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -msse2")

#subdirs
add_subdirectory(unittest)
