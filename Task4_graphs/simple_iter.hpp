#pragma once
#include "dense_CSR.hpp"

namespace simple_iter{
    using vector = std::vector<double>;

    std::pair<vector, unsigned> Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration); 

    std::pair<vector, unsigned> Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);

    std::pair<vector, unsigned> simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);
}