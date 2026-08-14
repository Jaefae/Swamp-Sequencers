# NOTE: no project()/cmake_minimum_required here — this is a subdirectory of the
# top-level SwampSequencer project and inherits its C++20 settings, exactly like
# src/. This directory is only added when SWAMP_BUILD_BENCH=ON, so the benchmark
# never costs anything in a default or test-only build.

add_executable(SwampBench "")
target_sources(SwampBench PRIVATE bench.cpp)

# SwampSeqLib carries the include directory as PUBLIC usage requirements.
target_link_libraries(SwampBench PRIVATE SwampSeqLib)

# -Wall/-Wextra/-Wpedantic are GCC/Clang spellings; MSVC needs /W4.
if(MSVC)
  target_compile_options(SwampBench PRIVATE /W4)
else()
  target_compile_options(SwampBench PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Deliberately no add_test(): the benchmark measures, it does not assert. CI runs
# it as an explicit step and uploads the CSV, so a slow machine can never turn
# into a red test run.
