#include <gtest/gtest.h>
#include "main.hpp"

#include <vector>
#include <string>
#include <cstring>

TEST(DRY_RUN, test_simple_prs_data)
{
    std::vector<std::string> args{
        "prs_tool",
        "--base", "data/Height.QC.gz",
        "--target", "data/EUR.QC",
        "--binary-target", "F",
        "--pheno", "data/EUR.height",
        "--cov", "data/EUR.covariate",
        "--base-maf", "MAF:0.01",
        "--base-info", "INFO:0.8",
        "--stat", "OR",
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
    int ret = main(argc, argv.data());

    // cleanup
    for (char* p : argv) {
        delete[] p;
    }

    EXPECT_EQ(ret, 0);
}
