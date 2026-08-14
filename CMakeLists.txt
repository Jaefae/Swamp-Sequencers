# Set minimum CMake version required
cmake_minimum_required(VERSION 3.15)

# Define project
project(SwampSequencer VERSION 1.0 LANGUAGES CXX)

# Define C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Centralize binary output for easier testing
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Build options 
# The TUI pulls in FTXUI and downloads the ~5 MB E. coli genome at configure
# time. CI only needs the library and the test suite, so both are opt-out:
#   cmake -B build -DSWAMP_BUILD_TUI=OFF
option(SWAMP_BUILD_TUI   "Build the FTXUI front-end and download the reference genome" ON)
option(SWAMP_BUILD_TESTS "Build the GoogleTest suite"                                  ON)
option(SWAMP_BUILD_BENCH "Build the suffix array vs. suffix tree benchmark"            OFF)
option(SWAMP_ENABLE_ASAN "Build with AddressSanitizer + UndefinedBehaviorSanitizer"    OFF)

# Sanitizers 
# Applied globally so the library, tests, and benchmark are all instrumented.
if(SWAMP_ENABLE_ASAN)
  if(MSVC)
    message(FATAL_ERROR "SWAMP_ENABLE_ASAN requires a GCC or Clang toolchain.")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address,undefined)
endif()

# Enable CTest
enable_testing()

# Fetch GTest framework
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.17.0
)

if(SWAMP_BUILD_TUI)
  FetchContent_Declare(ftxui
    GIT_REPOSITORY https://github.com/ArthurSonzogni/ftxui
    GIT_TAG v6.1.9
  )
  FetchContent_MakeAvailable(ftxui)
endif()

# Add src subdirectory
add_subdirectory(src)

if(SWAMP_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(SWAMP_BUILD_TESTS)
  # Option for Windows compatibility
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_subdirectory(tests)
endif()
