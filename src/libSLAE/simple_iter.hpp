#pragma once
#include "dense_CSR.hpp"
#include <cmath>
#define _USE_MATH_DEFINES

namespace simple_iter{
    using vector = std::vector<double>;

    vector Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration); 

    vector Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);

    vector simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);

    //Возвращает 2^n корней на отрезке [lambda_min, lambda_max] в оптимальном порядке
    vector find_Chebyshev_roots(unsigned n, double lambda_min, double lambda_max);
}