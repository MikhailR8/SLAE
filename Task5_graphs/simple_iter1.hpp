#pragma once
#include "dense_CSR.hpp"
#include <cmath>
#include <chrono>
#include <fstream>
#define _USE_MATH_DEFINES

namespace simple_iter{
    using steady_clock = std::chrono::steady_clock;
    using vector = std::vector<double>;
    class Timer {
    public:
        Timer(std::fstream& w) : start(steady_clock::now()), writer(w) {}
        ~Timer() {
            writer << std::chrono::duration_cast<std::chrono::microseconds>(steady_clock::now()
                - start).count()
                << std::endl;
        }
    private:
        steady_clock::time_point start;
        std::fstream& writer;
    };

    void Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration); 

    void Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);

    void simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration);

    void simple_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double lambda_min, double lambda_max,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration, unsigned quantity_roots);

    //Возвращает 2^n корней на отрезке [lambda_min, lambda_max] в оптимальном порядке
    vector find_Chebyshev_roots(unsigned n, double lambda_min, double lambda_max);
}