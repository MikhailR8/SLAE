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
        //Корни генерируются последовательно, корни парные - к каждому из первой половины найдётся пара
        //во второй с другим знаком, но порядок их неодинаков, так что для соблюдения порядка честно 
        //генерируем все
        for (auto i = 1u; i < count; i++){
            auto cos_now = root;
            root = root * cos_const - sin_now * sin_const;
            sin_now = sin_now * cos_const + sin_const * cos_now;
            roots[indices[i]] = root;
        }
        //Меняем местами корни и растягиваем, получится красиво: корни расставятся по строгому убыванию
        //(Честно говоря, можно было и sort() использовать, алгоритм выше тоже примерно за (n log n)
        //индексы расставляет, хотя там нет проверок, так что sort всё же помедленнее будет)
        vector res(count);
        for(auto i = 0u; i < count; i++){
            res[i] = 1 / ((lambda_max + lambda_min) / 2.0 +
             (lambda_max - lambda_min) / 2.0 * roots[indices[i]]);
        }
        return res;
    }
}