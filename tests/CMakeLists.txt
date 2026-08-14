add_executable(SwampTests)

target_sources(SwampTests PRIVATE
  swamp_test.cpp
)

target_link_libraries(SwampTests PRIVATE
  GTest::gtest_main
  SwampSeqLib
)

include(GoogleTest)
gtest_discover_tests(SwampTests)
