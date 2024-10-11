cmake_minimum_required(VERSION 3.5)
project(CvxCompress C CXX)

# Ensure C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for the 'simde' submodule
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/simde/x86/avx.h")
    message(FATAL_ERROR "The 'simde' submodule is missing or incomplete. Please run 'git submodule update --init --recursive'")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR})

# Find OpenMP for C++ only
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Set initial compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -g ${OpenMP_CXX_FLAGS}")

# Check if we're on x86_64 architecture and add -mavx if supported
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    # Check if compiler supports -mavx
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_MAVX)
    if(COMPILER_SUPPORTS_MAVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    endif()
endif()

# Create object library for Wavelet_Transform_Slow (used by CvxCompress_GenCode)
add_library(gen_code_objects OBJECT
    Wavelet_Transform_Slow.cpp
)

# Build CvxCompress_GenCode executable
add_executable(CvxCompress_GenCode
    CvxCompress_GenCode.cpp
    $<TARGET_OBJECTS:gen_code_objects>
)
target_link_libraries(CvxCompress_GenCode OpenMP::OpenMP_CXX m)

# Generate Ds79_Base.cpp and Us79_Base.cpp into the source directory
add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/Ds79_Base.cpp ${CMAKE_SOURCE_DIR}/Us79_Base.cpp
    COMMAND CvxCompress_GenCode
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS CvxCompress_GenCode
    COMMENT "Generating Ds79_Base.cpp and Us79_Base.cpp in source directory"
)

add_custom_target(generate_base_files
    DEPENDS ${CMAKE_SOURCE_DIR}/Ds79_Base.cpp ${CMAKE_SOURCE_DIR}/Us79_Base.cpp
)

# Create object library for Wavelet_Transform_Fast
add_library(Wavelet_Transform_Fast OBJECT
    Wavelet_Transform_Fast.cpp
)
add_dependencies(Wavelet_Transform_Fast generate_base_files)

# Create object library for common sources
add_library(common_objects OBJECT
    CvxCompress.cpp
    Run_Length_Encode_Slow.cpp
    Block_Copy.cpp
    Read_Raw_Volume.cpp
    Wavelet_Transform_Slow.cpp
)

# Create shared library 'cvxcompress'
add_library(cvxcompress SHARED
    $<TARGET_OBJECTS:common_objects>
    $<TARGET_OBJECTS:Wavelet_Transform_Fast>
)
target_link_libraries(cvxcompress OpenMP::OpenMP_CXX m)

# Install the cvxcompress library
install(TARGETS cvxcompress
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install specific header files
install(FILES
    ${CMAKE_SOURCE_DIR}/CvxCompress.hxx
    ${CMAKE_SOURCE_DIR}/Wavelet_Transform_Slow.hxx
    ${CMAKE_SOURCE_DIR}/Run_Length_Encode_Slow.hxx
    ${CMAKE_SOURCE_DIR}/Run_Length_Escape_Codes.hxx
    ${CMAKE_SOURCE_DIR}/Read_Raw_Volume.hxx
    ${CMAKE_SOURCE_DIR}/Wavelet_Transform_Fast.hxx
    ${CMAKE_SOURCE_DIR}/Block_Copy.hxx
    DESTINATION include/cvxcompress
)

# Executable: CvxCompress_Test (static linking)
add_executable(CvxCompress_Test
    CvxCompress_Test.cpp
    $<TARGET_OBJECTS:common_objects>
    $<TARGET_OBJECTS:Wavelet_Transform_Fast>
)
target_link_libraries(CvxCompress_Test OpenMP::OpenMP_CXX m)

# Executable: CvxCompress_Test_Dyn (dynamic linking)
add_executable(CvxCompress_Test_Dyn
    CvxCompress_Test.cpp
)
add_dependencies(CvxCompress_Test_Dyn cvxcompress)
target_link_libraries(CvxCompress_Test_Dyn cvxcompress OpenMP::OpenMP_CXX m)

# Executable: Test_Compression
add_executable(Test_Compression
    Test_Compression.cpp
    $<TARGET_OBJECTS:common_objects>
    $<TARGET_OBJECTS:Wavelet_Transform_Fast>
)
target_link_libraries(Test_Compression OpenMP::OpenMP_CXX m)

# Executable: Compress_SEAM_Basin
add_executable(Compress_SEAM_Basin
    Compress_SEAM_Basin.cpp
)
add_dependencies(Compress_SEAM_Basin cvxcompress)
target_link_libraries(Compress_SEAM_Basin cvxcompress OpenMP::OpenMP_CXX m)

# Executable: Test_With_Generated_Input
add_executable(Test_With_Generated_Input
    Test_With_Generated_Input.cpp
)
add_dependencies(Test_With_Generated_Input cvxcompress)
target_link_libraries(Test_With_Generated_Input cvxcompress OpenMP::OpenMP_CXX m)