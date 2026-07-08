include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(nurbslab_tests 
    test_evaluate.cpp
    test_projection.cpp
    test_intersection.cpp
)

target_link_libraries(nurbslab_tests PRIVATE nurbslab GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(nurbslab_tests)