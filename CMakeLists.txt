cmake_minimum_required(VERSION 3.10)

project(DarkChess)

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the include directory to the project
include_directories(include)

file(GLOB C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/KZG_BN254/*.c")

add_library(c_lib STATIC ${C_SOURCES})

target_include_directories(c_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/KZG_BN254
)

# Find Botan
find_package(PkgConfig REQUIRED)
pkg_check_modules(BOTAN REQUIRED botan-3)

# Explicitly set Botan library path
set(BOTAN_LIBRARY_DIR "/opt/homebrew/Cellar/botan/3.5.0/lib")

# Add the executable
add_executable(DarkChess src/main.cpp src/utils.cpp src/Chess.cpp src/BiVOT.cpp)

# Specify the Icicle directories
set(ICICLE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/icicle/icicle/include)
set(ICICLE_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/icicle/build)

# Find the Icicle libraries
find_library(ICICLE_FIELD_LIBRARY
    NAMES icicle_field_bn254
    PATHS ${ICICLE_LIB_DIR}
    NO_DEFAULT_PATH
)

find_library(ICICLE_DEVICE_LIBRARY
    NAMES icicle_device
    PATHS ${ICICLE_LIB_DIR}
    NO_DEFAULT_PATH
)

find_library(ICICLE_CURVE_LIBRARY
    NAMES icicle_curve_bn254
    PATHS ${ICICLE_LIB_DIR}
    NO_DEFAULT_PATH
)

find_library(ICICLE_HASH_LIBRARY
    NAMES icicle_hash
    PATHS ${ICICLE_LIB_DIR}
    NO_DEFAULT_PATH
)

# Check if the libraries were found
if(NOT ICICLE_FIELD_LIBRARY OR NOT ICICLE_DEVICE_LIBRARY OR NOT ICICLE_CURVE_LIBRARY OR NOT ICICLE_HASH_LIBRARY)
    message(FATAL_ERROR "Icicle libraries not found. Please ensure the Icicle library is built.")
endif()

# Include Icicle headers
target_include_directories(DarkChess PRIVATE ${ICICLE_INCLUDE_DIR})

# Include Botan headers
target_include_directories(DarkChess PRIVATE ${BOTAN_INCLUDE_DIRS})

# Link against Icicle libraries and Botan
target_link_libraries(DarkChess PRIVATE 
    ${ICICLE_FIELD_LIBRARY}
    ${ICICLE_DEVICE_LIBRARY}
    ${ICICLE_CURVE_LIBRARY}
    ${ICICLE_HASH_LIBRARY}
    c_lib
)


# Explicitly link against Botan
target_link_libraries(DarkChess PRIVATE ${BOTAN_LIBRARY_DIR}/libbotan-3.dylib)

# Add Botan compile options
target_compile_options(DarkChess PRIVATE ${BOTAN_CFLAGS_OTHER})

# Add additional compiler flags for Botan compatibility
target_compile_options(DarkChess PRIVATE -DBOTAN_ENABLE_STD_RANGES -DBOTAN_FORCE_STD_CONCEPTS)

# Print Botan information for debugging
message(STATUS "Botan include dirs: ${BOTAN_INCLUDE_DIRS}")
message(STATUS "Botan libraries: ${BOTAN_LIBRARIES}")
message(STATUS "Botan cflags: ${BOTAN_CFLAGS_OTHER}")
message(STATUS "Botan library dir: ${BOTAN_LIBRARY_DIR}")

# Print compiler information for debugging
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")