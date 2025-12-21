include(FetchContent)

# Invoke Gtest
fetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Set up a run from the tutorial: https://choishingwan.github.io/PRS-Tutorial/
add_executable(test_tutorial_run test_tutorial_run.cpp)
target_link_libraries(test_tutorial_run GTest::gtest_main 
                                        bgen 
                                        gzstream 
                                        plink 
                                        prsice_lib 
                                        genotyping 
                                        regression 
                                        utility coverage_config)
target_include_directories(test_tutorial_run PUBLIC ${CMAKE_SOURCE_DIR}/inc)
add_test(NAME test_tutorial_run COMMAND test_tutorial_run)