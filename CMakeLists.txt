cmake_minimum_required(VERSION 3.10)
project(DarkChess)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the include directory to the project
include_directories(include)

# Add the executable
add_executable(DarkChess src/main.cpp src/utils.cpp)

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

# Link against Icicle libraries
target_link_libraries(DarkChess PRIVATE 
    ${ICICLE_FIELD_LIBRARY}
    ${ICICLE_DEVICE_LIBRARY}
    ${ICICLE_CURVE_LIBRARY}
    ${ICICLE_HASH_LIBRARY}
)