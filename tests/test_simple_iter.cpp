#include <gtest/gtest.h>
#include <libSLAE/simple_iter.hpp>

TEST(Matrix_test, first_matrix) {
    std::vector<std::vector<double>> arr = 
    std::vector<std::vector<double>>({{6, -2, 2}, {-3, 6, -1}, {-2, 5, 8}});
    auto b = std::vector<double>({5, 0, -4});
    auto x_0 = std::vector<double>({-4, 4, 4});
    auto mat = dense_CSR::Matrix_CSR(arr);
    auto out = simple_iter::Jacobi_iter(mat, x_0, b, 0.01, 1, 20);
    auto out1 = simple_iter::Gauss_Seidel_iter(mat, x_0, b, 0.01, 1, 20);
    //tau_opt = 2 / 12.214 = 0.1637
    auto out2 = simple_iter::simple_iter(mat, x_0, b, 0.1637, 0.01, 1, 20); 
    ASSERT_NEAR(0.50, out[1], 0.01);
    ASSERT_NEAR(0.50, out1[1], 0.01);
    ASSERT_NEAR(0.50, out2[1], 0.01);
}

int main(int argc, char **argv) {
// int main(){
//     // https://www.geogebra.org/m/XBKbmXY7
//     std::vector<std::vector<double>> arr = 
//     std::vector<std::vector<double>>({{6, -2, 2}, {-3, 6, -1}, {-2, 5, 8}});
//     auto b = std::vector<double>({5, 0, -4});
//     auto x_0 = std::vector<double>({-4, 4, 4});
//     auto mat = dense_CSR::Matrix_CSR(arr);
//     auto out = simple_iter::Jacobi_iter(mat, x_0, b, 0.01, 1, 20);
//     auto out1 = simple_iter::Gauss_Seidel_iter(mat, x_0, b, 0.01, 1, 20);
//     //tau_opt = 2 / 12.214 = 0.1637
//     auto out2 = simple_iter::simple_iter(mat, x_0, b, 0.1637, 0.01, 1, 20); 
//     dense_CSR::print_vector(out);
//     dense_CSR::print_vector(out1);
//     dense_CSR::print_vector(out2);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

