#include "../src/libSLAE/Iter_solvers.hpp"
#include "../src/libSLAE/GMRES.hpp"
#include <random>


int main(){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist2(1,10);

    //считываем матрицу, сгенерированную three_diag_generator.py
    unsigned n = 100u;
    std::vector<vector> data(n, vector(n));
    vector vec(n);
    vector b(n);
    std::ifstream in("../matrix.txt");
    double lambda_min, lambda_max;
    if (in.is_open()) {
        for(auto i = 0u; i < n; i++){
            for(auto j = 0u; j < n; j++){
                in >> data[i][j];
            }
        }
        in >> lambda_min >> lambda_max;
    }
    in.close();

    //Заполняем вектор
    for (auto i = 0u; i < n; i++){
        vec[i] = static_cast<double>(dist2(rng));
    }


    for (auto i = 0u; i < n; i++){
        b[i] = static_cast<double>(dist2(rng));
    }

    double tau_opt = 2 / (lambda_min + lambda_max);
    auto matrix = dense_CSR::Matrix_CSR(data);
    double target_discrepancy = 10E-9;
    unsigned max_iteration = 100u;
    // auto vec1 = Iter_solvers::Jacobi_iter(matrix, vec, b, target_discrepancy, 1, max_iteration, true);
    // auto vec2 = Iter_solvers::Gauss_Seidel_iter(matrix, vec, b, target_discrepancy, 1, max_iteration, true);    
    // auto vec3 = Iter_solvers::Symmetric_Gauss_Seidel_iter(matrix, vec, b, target_discrepancy,
    //  1, max_iteration, true);
    // auto vec4 = Iter_solvers::Symmetric_Gauss_Seidel_iter_boost(matrix, vec, b, target_discrepancy,
    //  1, max_iteration, 0.0001, true);
    // auto vec5 = Iter_solvers::simple_iter(matrix, vec, b, tau_opt, target_discrepancy, 1,
    //  max_iteration, true);
    // auto vec6 = Iter_solvers::simple_iter_boost(matrix, vec, b, lambda_min, lambda_max,
    //  target_discrepancy, 1, max_iteration, 4, true);  
    // auto vec7 = Iter_solvers::fastest_descent(matrix, vec, b, target_discrepancy, 1, max_iteration, true);
    auto vec8 = Iter_solvers::Conjugate_gradient(matrix, vec, b,
     target_discrepancy, max_iteration, 10E-20, true); 
    auto vec9 = Iter_solvers::GMRES(matrix, vec, b, target_discrepancy, max_iteration, true); 
    dense_CSR::print_vector(vec8); 
    std::cout << std::endl;
    dense_CSR::print_vector(vec9); 
    std::cout << std::endl;
    // dense_CSR::print_vector(vec3); 
    // std::cout << std::endl;
    // dense_CSR::print_vector(vec4); 
    // std::cout << std::endl;
    // dense_CSR::print_vector(vec5); 
    // std::cout << std::endl;
    // dense_CSR::print_vector(vec6); 
    // std::cout << std::endl;
    // dense_CSR::print_vector(vec7); 
    // std::cout << std::endl;
    // dense_CSR::print_vector(vec7); 
}