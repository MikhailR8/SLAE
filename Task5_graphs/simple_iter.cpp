#include "simple_iter.hpp"

namespace simple_iter{
    using dense_CSR::operator+;
    using dense_CSR::operator*;
    using dense_CSR::operator-;

    vector Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy=0.1, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        vector discrepancy(max_iteration);
        for (auto i = 0u; i < max_iteration; i++){
            discrepancy[i] = dense_CSR::get_length(A * res - b);
            for(auto j = 0u; j < 10u; j++){
                res = A.diag_reverse_multiplication(b - A.LU_multiplication(res));
            }
        }
        return discrepancy;
    }

    vector Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy=0.1, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        vector discrepancy(max_iteration);
        for (auto i = 0u; i < max_iteration; i++){
            discrepancy[i] = dense_CSR::get_length(A * res - b);
            for(auto j = 0u; j < 10u; j++){
                res = A.Gauss_Seidel_iteration(res, b);
            }
        }
        return discrepancy;
    }

    vector simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy=0.1, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        vector discrepancy(max_iteration);
        for (auto i = 0u; i < max_iteration; i++){
            discrepancy[i] = dense_CSR::get_length(A * res - b);
            for(auto j = 0u; j < 10u; j++){
                res = res - tau * (A * res - b);
            }
        }
        return discrepancy;
    }   

    vector simple_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double lambda_min, double lambda_max, double target_discrepancy=0.1,
    unsigned frequency_checking=3u, unsigned max_iteration=100u, unsigned quantity_roots = 3u){
        auto res = vector(x_0);
        auto steps = find_Chebyshev_roots(quantity_roots, lambda_min, lambda_max);
        auto step = 0u;
        unsigned count = static_cast<unsigned>(std::pow(2, quantity_roots));
        vector discrepancy(max_iteration);
        for (auto i = 0u; i < max_iteration; i++){
            discrepancy[i] = dense_CSR::get_length(A * res - b);
            for(auto j = 0u; j < 10u; j++){
                step++;
                res = res - steps[step % count] * (A * res - b);
            }
        }
        return discrepancy;
    }  

    vector find_Chebyshev_roots(unsigned n, double lambda_min, double lambda_max){
        unsigned count = static_cast<unsigned>(std::pow(2, n));
        std::vector<unsigned> indices(count);
        vector roots(count);

        //Заполняем массив индексов корней
        indices[count / 2] = 1;
        auto step = count / 2u;
        for(auto degree = 2u; degree <= n; degree++){
            step /= 2u;
            auto parameter = static_cast<unsigned>(std::pow(2, degree));
            for(auto place = step; place < count; place += (2 * step)){
                indices[place] = parameter - 1 - indices[place - step]; 
            }
        }
        
        double root = std::cos(M_PI / 2.0 / std::pow(2, n));
        double sin_now = std::sin(M_PI / 2.0 / std::pow(2, n));
        double sin_const = 2 * root * sin_now;
        double cos_const = (root * root) - (sin_now * sin_now);
        roots[0] = root;
        roots[count - 1u] = -1.0 * root;

        for (auto i = 1u; i < (count / 2u); i++){
            auto cos_now = root;
            root = root * cos_const - sin_now * sin_const;
            sin_now = sin_now * cos_const + sin_const * cos_now;
            roots[i] = root;
            roots[count - 1u - i] = -1.0 * root;
        }

        vector res(count);
        for(auto i = 0u; i < count; i++){
            res[i] = 1.0 / ((lambda_max + lambda_min) / 2.0
             + (lambda_max - lambda_min) / 2.0 *  roots[indices[i]]);
        }

        return res;
    }
}