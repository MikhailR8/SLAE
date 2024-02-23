#include <gtest/gtest.h>
#include <libSLAE/Householder_alg.hpp>
#include <iostream>
#include <vector>

TEST(Matrix_test, first_matrix) {
    auto mat = dense_CSR::Matrix({{3, 5, 8}, {5, 3, 7}, {2, 4, 9}});
    auto solution = Householder::solve_qr(std::vector<double>({8, 7, 6}), mat);
    ASSERT_NEAR(solution[0], 0.717, 0.001);
    ASSERT_NEAR(solution[1], 1.239, 0.001);
    ASSERT_NEAR(solution[2], -0.043, 0.001);
}

int main(int argc, char **argv){
    // auto mat = dense_CSR::Matrix({{3, 5, 8}, {5, 3, 7}, {2, 4, 9}});
    // auto pair = Householder::QR(mat);
    // auto Q_ref = dense_CSR::Matrix({{0.49, 0.58, -0.66}, {0.81, -0.58, 0.09}, {0.32, 0.58, 0.75}});
    // auto R_ref = dense_CSR::Matrix({{6.16, 6.16, 12.49}, {0, 3.46, 5.77}, {0, 0, 2.15}});
    // dense_CSR::print_matrix<dense_CSR::Matrix>(pair.first);
    // std::cout << " " << std::endl;
    // dense_CSR::print_matrix<dense_CSR::Matrix>(pair.second);
    // std::cout << "Refs:" << std::endl;
    // dense_CSR::print_matrix<dense_CSR::Matrix>(Q_ref);
    // std::cout << " " << std::endl;
    // dense_CSR::print_matrix<dense_CSR::Matrix>(R_ref);
    // auto solution = Householder::solve_qr(std::vector<double>({8, 7, 6}), mat);
    // std::cout << solution[0] << " " << solution[1] << " " << solution[2] << std::endl;
    // std::cout << "Refs:" << std::endl;
    // std::cout << 0.717 << " " << 1.239 << " " << -0.043 << std::endl;
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}