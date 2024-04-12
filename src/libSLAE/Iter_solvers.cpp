#include "Iter_solvers.hpp"

namespace Iter_solvers{
    using dense_CSR::operator+;
    using dense_CSR::operator*;
    using dense_CSR::operator-;

    vector Jacobi_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});

        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy = dense_CSR::get_length(A * res - b);
            if(discrepancy < target_discrepancy){
                if(testmode) print_to_file(out, "Jacobi");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy);
                        {
                            Timer T(&out.second);
                            res = A.Jacobi_iteration(res, b);
                        }
                    }
                    else{
                        res = A.Jacobi_iteration(res, b);
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Jacobi");
        return res;
    }

    vector Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});

        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy = dense_CSR::get_length(A * res - b);
            if(discrepancy < target_discrepancy){
                if(testmode) print_to_file(out, "Gauss_Seidel");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy);
                        {
                            Timer T(&out.second);
                            res = A.Gauss_Seidel_iteration(res, b);
                        }
                    }
                    else{
                        res = A.Gauss_Seidel_iteration(res, b);
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Gauss_Seidel");
        return res;
    }

    vector Symmetric_Gauss_Seidel_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration,
    bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});

        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy = dense_CSR::get_length(A * res - b);
            if(discrepancy < target_discrepancy){
                if(testmode) print_to_file(out, "Symmetric_Gauss_Seidel");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy);
                        {
                            Timer T(&out.second);
                            res = A.Symmetric_Gauss_Seidel_iteration(res, b);
                        }
                    }
                    else{
                        res = A.Symmetric_Gauss_Seidel_iteration(res, b);
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Symmetric_Gauss_Seidel");
        return res;
    }

    vector Symmetric_Gauss_Seidel_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double target_discrepancy, unsigned frequency_checking,
    unsigned max_iteration, double rho, bool testmode){
        test_pair out({vector(), std::vector<unsigned long long>()});

        double mu_i_minus_1 = 1.0; //mu_0
        double mu_i = 1.0 / rho; //mu_1
        double mu_i_plus_1 = 2.0 / rho * mu_i - mu_i_minus_1; //mu_2
        vector y_i(x_0); //y_0
        vector y_i_plus_1(b.size()); //y_1
        if(testmode) {
            {
                Timer T(&out.second);
                y_i_plus_1 = A.Symmetric_Gauss_Seidel_iteration(y_i, b);
            }
        } else{y_i_plus_1 = A.Symmetric_Gauss_Seidel_iteration(y_i, b);}

    for (auto i = 0u; i < max_iteration; i+=frequency_checking){
        double discrepancy = dense_CSR::get_length(A * y_i_plus_1 - b);
        if(discrepancy < target_discrepancy){
            if(testmode) print_to_file(out, "Symmetric_Gauss_Seidel_boost");
            return y_i_plus_1;
        }
        else{
            for(auto j = 0u; j < frequency_checking; j++){
                if(testmode){
                out.first.push_back(discrepancy);
                {
                    Timer T(&out.second);
                    auto temp_y = y_i_plus_1;
                    y_i_plus_1 = (2.0 * mu_i / rho / mu_i_plus_1) * 
                    (A.Symmetric_Gauss_Seidel_iteration(y_i_plus_1, b))
                    - (mu_i_minus_1 / mu_i_plus_1) * (y_i);
                    y_i = temp_y;

                    auto temp_mu_1 = mu_i;
                    auto temp_mu_2 = mu_i_plus_1;
                    mu_i_minus_1 = mu_i;
                    mu_i = mu_i_plus_1;
                    mu_i_plus_1 = 2.0 / rho * temp_mu_2 - temp_mu_1;
                }
            }
                else{
                auto temp_y = y_i_plus_1;
                y_i_plus_1 = (2.0 * mu_i / rho / mu_i_plus_1) * 
                (A.Symmetric_Gauss_Seidel_iteration(y_i_plus_1, b))
                 - (mu_i_minus_1 / mu_i_plus_1) * (y_i);
                y_i = temp_y;

                auto temp_mu_1 = mu_i;
                auto temp_mu_2 = mu_i_plus_1;
                mu_i_minus_1 = mu_i;
                mu_i = mu_i_plus_1;
                mu_i_plus_1 = 2.0 / rho * temp_mu_2 - temp_mu_1;
            }
            }
        }
    }
    if(testmode) print_to_file(out, "Symmetric_Gauss_Seidel_boost");
    return y_i_plus_1;     
    }   

    vector simple_iter(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b, double tau,
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration, bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});

        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy = dense_CSR::get_length(A * res - b);
            if(discrepancy < target_discrepancy){
                if(testmode) print_to_file(out, "Simple_iter");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy);
                        {
                            Timer T(&out.second);
                            res = res - tau * (A * res - b);
                        }
                    }
                    else{
                        res = res - tau * (A * res - b);
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Simple_iter");
        return res;
    } 
    
    vector simple_iter_boost(const dense_CSR::Matrix_CSR& A, const vector& x_0,
    const vector& b, double lambda_min, double lambda_max, double target_discrepancy,
    unsigned frequency_checking, unsigned max_iteration, unsigned quantity_roots, bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});
        vector steps(1u << quantity_roots);

        if(testmode){
            {
                Timer(&out.second);
                steps = find_Chebyshev_roots(quantity_roots, lambda_min, lambda_max);
            }
        } else{steps = find_Chebyshev_roots(quantity_roots, lambda_min, lambda_max);}
        auto step = 0u;
        unsigned count = 1u << quantity_roots;
        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy = dense_CSR::get_length(A * res - b);
            if(discrepancy < target_discrepancy){
                if(testmode) print_to_file(out, "Simple_iter_boost");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy);
                        {
                            Timer T(&out.second);
                            res = res - steps[step % count] * (A * res - b);
                            step++;
                        }
                    }
                    else{
                        res = res - steps[step % count] * (A * res - b);
                        step++;
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Simple_iter_boost");
        return res;


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
    double target_discrepancy, unsigned frequency_checking, unsigned max_iteration, bool testmode){
        auto res = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});
        vector discrepancy = A * res - b;

        for (auto i = 0u; i < max_iteration; i+=frequency_checking){
            double discrepancy_mod = dense_CSR::get_length(A * res - b);
            if(discrepancy_mod < target_discrepancy){
                if(testmode) print_to_file(out, "Fastest_descent");
                return res;
            }
            else{
                for(auto j = 0u; j < frequency_checking; j++){
                    if(testmode){
                        out.first.push_back(discrepancy_mod);
                        {
                            Timer T(&out.second);
                            double tau = (discrepancy * discrepancy) / (discrepancy * (A * discrepancy));
                            res = res - tau * (discrepancy);
                            discrepancy = A * res - b;
                        }
                    }
                    else{
                        double tau = (discrepancy * discrepancy) / (discrepancy * (A * discrepancy));
                        res = res - tau * (discrepancy);
                        discrepancy = A * res - b;
                    }
                }
            }
        }
        if(testmode) print_to_file(out, "Fastest_descent");
        return res;
    }

    vector Conjugate_gradient(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned max_iteration, double epsilon, bool testmode){
        //Матрица A O(n*) memory
        auto res = vector(x_0); //x_0 O(n) memory
        test_pair out({vector(), std::vector<unsigned long long>()});
        vector discrepancy = A * res - b; //r_0 O(n) memory
        vector d = discrepancy; //d_0 O(n) memory
        double discrepancy_mod = dense_CSR::get_length(discrepancy); //O(const) memory
        double d_mod = discrepancy_mod; //O(n) memory
        double alpha = 0.0; //O(const) memory
        double beta = 0.0; //O(const) memory
        vector discrepancy_plus_1(x_0.size()); //O(n) memory

        if(testmode) out.first.push_back(discrepancy_mod);
        if(discrepancy_mod < target_discrepancy){
            if(testmode) print_to_file(out, "Conjugate_gradient");
            return res;
        }

        for (auto i = 0u; i < max_iteration; i++){
            if(testmode){
                {
                    Timer T(&out.second);
                    alpha = (discrepancy * discrepancy) / (d * (A * d));
                    res = res - alpha * d;
                    discrepancy_plus_1 = A * res - b;
                    beta = (discrepancy_plus_1 * discrepancy_plus_1) / (discrepancy * discrepancy);
                    d = discrepancy_plus_1 + beta * d;
                    discrepancy = discrepancy_plus_1;
                    discrepancy_mod = dense_CSR::get_length(discrepancy_plus_1);
                    d_mod = dense_CSR::get_length(d);           
                }
                out.first.push_back(discrepancy_mod);
                if(discrepancy_mod < target_discrepancy || d_mod < epsilon){
                    print_to_file(out, "Conjugate_gradient");
                    return res;                   
                }                            
            }
            else{
                alpha = (discrepancy * discrepancy) / (d * (A * d)); //O(n*) speed
                res = res - alpha * d; //O(n) speed
                discrepancy_plus_1 = A * res - b; //O(n*) speed
                beta = (discrepancy_plus_1 * discrepancy_plus_1) / (discrepancy * discrepancy); //O(n) speed
                d = discrepancy_plus_1 + beta * d; //O(n) speed
                discrepancy = discrepancy_plus_1; //O(n) speed
                discrepancy_mod = dense_CSR::get_length(discrepancy_plus_1); //O(n) speed
                d_mod = dense_CSR::get_length(d); //O(n) speed

                // dense_CSR::print_vector(discrepancy);                
                // dense_CSR::print_vector(d);
                
                // std::cout << std::endl;
                if(discrepancy_mod < target_discrepancy || d_mod < epsilon){
                    return res;                   
                }             
            }
        }
        if(testmode) print_to_file(out, "Conjugate_gradient");
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