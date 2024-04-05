#pragma once
#include "dense_CSR.hpp"
#include "Test_speed_tools.hpp"
#include <cmath>
#define _USE_MATH_DEFINES

namespace Iter_solvers{
    using vector = std::vector<double>;

    vector Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode = false); 

    vector Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode = false);

    vector Symmetric_Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode = false);

    vector Symmetric_Gauss_Seidel_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double target_discrepancy, unsigned frequency_checking,
    unsigned max_iteration, double rho = 0.8, bool testmode = false);

    vector simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode = false);

    vector simple_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double lambda_min, double lambda_max, double target_discrepancy,
    unsigned frequency_checking, unsigned max_iteration, unsigned quantity_roots, bool testmode = false);

    vector fastest_descent(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode = false);

    //Возвращает 2^n корней на отрезке [lambda_min, lambda_max] в оптимальном порядке
    vector find_Chebyshev_roots(unsigned n, double lambda_min, double lambda_max);
}