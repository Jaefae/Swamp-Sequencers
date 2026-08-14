# Only added when SWAMP_BUILD_BENCH=ON; inherits the top-level C++20 settings.

add_executable(SwampBench "")
target_sources(SwampBench PRIVATE bench.cpp)

target_link_libraries(SwampBench PRIVATE SwampSeqLib)

if(MSVC)
  target_compile_options(SwampBench PRIVATE /W4)
else()
  target_compile_options(SwampBench PRIVATE -Wall -Wextra -Wpedantic)
endif()

# No add_test(): the benchmark measures, it does not assert, so a slow machine
# must not turn into a failing test run.
