#include "simple_iter.hpp"

namespace simple_iter{
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
                    res = A.diag_reverse_multiplication(b - A.LU_multiplication(res));
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
}