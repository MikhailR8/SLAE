#include "Iter_solvers.hpp"

namespace Iter_solvers{
    using dense_CSR::operator+;
    using dense_CSR::operator*;
    using dense_CSR::operator-;

    vector Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy=5.0, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(A * res - b) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    res = A.Jacobi_iteration(res, b);
                }
            }
        }
        return res;
    }

    vector Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy=5.0, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(A * res - b) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    res = A.Gauss_Seidel_iteration(res, b);
                }
            }
        }
        return res;
    }

    vector Symmetric_Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double target_discrepancy=5.0,
    unsigned frequency_checking=5u, unsigned max_iteration=100u) {
        auto res = vector(x_0);
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(A * res - b) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    res = A.Symmetric_Gauss_Seidel_iteration(res, b);
                }
            }
        }
        return res;        
    }

    vector simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy=5.0, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(A * res - b) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    res = res - tau * (A * res - b);
                }
            }
        }
        return res;
    } 
    
    vector simple_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double lambda_min, double lambda_max, double target_discrepancy=5.0,
    unsigned frequency_checking=5u, unsigned max_iteration=100u, unsigned quantity_roots = 3u){
        auto res = vector(x_0);
        auto steps = find_Chebyshev_roots(quantity_roots, lambda_min, lambda_max);
        auto step = 0u;
        unsigned count = 1u << quantity_roots;
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(A * res - b) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    res = res - steps[step % count] * (A * res - b);
                    step++;
                }
            }
        }
        return res;
    }    

    vector fastest_descent(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy=5.0, unsigned frequency_checking=5u, unsigned max_iteration=100u){
        auto res = vector(x_0);
        vector discrepancy = A * res - b;
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            if(dense_CSR::get_length(discrepancy) < target_discrepancy){
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    double tau = (discrepancy * discrepancy) / (discrepancy * (A * discrepancy));
                    res = res - tau * (discrepancy);
                    discrepancy = A * res - b;
                }
            }
        }
        return res;
    }

    vector find_Chebyshev_roots(unsigned n, double lambda_min, double lambda_max){
        unsigned count = 1u << n; //2^n
        std::vector<unsigned> indices(count);
        vector roots(count);

        //Заполняем массив индексов корней
        indices[count / 2] = 1;
        auto step = count / 2u;
        for(auto degree = 2u; degree <= n; degree++){
            step /= 2u;
            unsigned parameter = 1u << degree;
            for(auto place = step; place < count; place += (2 * step)){
                indices[place] = parameter - 1 - indices[place - step]; 
            }
        }
        
        double root = std::cos(M_PI / 2.0 / (1u << n));
        double sin_now = std::sin(M_PI / 2.0 / (1u << n));
        double sin_const = 2 * root * sin_now;
        double cos_const = (root * root) - (sin_now * sin_now);
        roots[0] = root;
        roots[count - 1u] = -1 * root;

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