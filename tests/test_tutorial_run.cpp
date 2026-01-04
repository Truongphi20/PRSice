#include <gtest/gtest.h>
#include "prs_score.hpp"

#include <vector>
#include <string>
#include <cstring>

#ifndef TEST_DATA_FOLDER
#define TEST_DATA_FOLDER
#endif

TEST(DRY_RUN, test_simple_prs_data)
{
    std::vector<std::string> args{
        "prs_tool",
        "--base", std::string(TEST_DATA_FOLDER) + "/Height.QC.gz",
        "--target", std::string(TEST_DATA_FOLDER) + "/EUR.QC",
        "--binary-target", "F",
        "--pheno", std::string(TEST_DATA_FOLDER) + "/EUR.height",
        "--cov", std::string(TEST_DATA_FOLDER) + "/EUR.covariate",
        "--base-maf", "MAF:0.01",
        "--base-info", "INFO:0.8",
        "--stat", "OR",
        "--thread", "1",
        "--or",
        "--out", "EUR"
    };

    std::vector<char*> argv;
    argv.reserve(args.size());

    for (const auto& s : args) {
        char* buf = new char[s.size() + 1];
        std::strcpy(buf, s.c_str());
        argv.push_back(buf);
    }

    int argc = static_cast<int>(argv.size());
    int ret = run_prs(argc, argv.data());
    
    EXPECT_EQ(ret, 0);
}
