#include <gtest/gtest.h>
#include <libSLAE/Iter_solvers.hpp>

// TEST(Matrix_test, first_matrix) {
//     std::vector<std::vector<double>> arr = 
//     std::vector<std::vector<double>>({{6, -2, 2}, {-3, 6, -1}, {-2, 5, 8}});
//     auto b = std::vector<double>({5, 0, -4});
//     auto x_0 = std::vector<double>({-4, 4, 4});
//     auto mat = dense_CSR::Matrix_CSR(arr);
//     auto out = Iter_solvers::Symmetric_Gauss_Seidel_iter(mat, x_0, b, 0.01, 1, 20);
//     auto out1 = Iter_solvers::fastest_descent(mat, x_0, b, 0.01, 1, 20); 
//     auto out2 = Iter_solvers::simple_iter_boost(mat, x_0, b, 3.2, 9, 0.01, 1, 20, 2);
//     auto out3 = Iter_solvers::Symmetric_Gauss_Seidel_iter_boost(mat, x_0, b, 0.01, 1, 20);
//     ASSERT_NEAR(0.50, out[1], 0.01);
//     ASSERT_NEAR(0.50, out1[1], 0.01);
//     ASSERT_NEAR(0.50, out2[1], 0.01);
// }

// int main(int argc, char **argv) {
int main(){
    // https://www.geogebra.org/m/XBKbmXY7
    std::vector<std::vector<double>> arr = 
    std::vector<std::vector<double>>({{10, -1, 2}, {-1, 20, -3}, {2, -3, 11}});
    auto b = std::vector<double>({5, 0, -4});
    auto x_0 = std::vector<double>({-4, 4, 4});
    auto mat = dense_CSR::Matrix_CSR(arr);
    // auto out = Iter_solvers::Jacobi_iter(mat, x_0, b, 0.01, 1, 20);
    // auto out1 = Iter_solvers::Symmetric_Gauss_Seidel_iter(mat, x_0, b, 0.01, 1, 20);
    // auto out2 = Iter_solvers::fastest_descent(mat, x_0, b, 0.01, 1, 20); 
    auto out3 = Iter_solvers::Symmetric_Gauss_Seidel_iter_boost(mat, x_0, b, 10E-200, 1, 2000, 0.5);
    auto out4 = Iter_solvers::Conjugate_gradient(mat, x_0, b, 10E-100, 2000, 10E-25);
    // dense_CSR::print_vector(out);
    // dense_CSR::print_vector(out1);
    // dense_CSR::print_vector(out2);
    dense_CSR::print_vector(out3);
    dense_CSR::print_vector(out4); 
    // testing::InitGoogleTest(&argc, argv);
    // return RUN_ALL_TESTS();
}

