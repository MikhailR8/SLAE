#include "dense_CSR.hpp"
#include <gtest/gtest.h>

TEST(Matrix_test, first_matrix) {
    std::vector<std::vector<double>> arr = 
    std::vector<std::vector<double>>({{0, 2, 0}, {2, 0, 4}, {1, 2, 0}});
    auto mat = cmake_example::Matrix(arr);
    auto mat_csr = cmake_example::Matrix_CSR(arr);
    std::vector<double> v = std::vector<double>({1, 5, 3});
    auto out1 = mat * v;
    auto out2 = mat_csr * v;
    ASSERT_NEAR(10.0, out1[0], 0.0001);
    ASSERT_NEAR(14.0, out1[1], 0.0001);
    ASSERT_NEAR(11.0, out1[2], 0.0001);
    ASSERT_NEAR(10.0, out2[0], 0.0001);
    ASSERT_NEAR(14.0, out2[1], 0.0001);
    ASSERT_NEAR(11.0, out2[2], 0.0001);
}

int main(int argc, char **argv) {
    // std::vector<std::vector<double>> arr = 
    // std::vector<std::vector<double>>({{0, 2, 0}, {2, 0, 4}, {0, 0, 0}});
    // auto mat = cmake_example::Matrix(arr);
    // auto mat_csr = cmake_example::Matrix_CSR(arr);
    // cmake_example::print_matrix<cmake_example::Matrix>(mat);
    // cmake_example::print_matrix<cmake_example::Matrix_CSR>(mat_csr);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

