#include "simple_iter1.hpp"
#include <chrono>
#include <fstream>
#include <random>
#include <iostream>

using vector = std::vector<double>;

int main(){

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist2(1,10);

    //считываем матрицу, сгенерированную three_diag_generator.py
    unsigned n = 100u;
    std::vector<vector> data(n, vector(n));
    vector vec(n);
    vector b(n);
    std::ifstream in("matrix.txt");
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
    // std::cout << lambda_min << " " << lambda_max << std::endl;

    // for(auto i = 0u; i < 100u; i++){
    //     for(auto j = 0u; j < 100u; j++){
    //             std::cout << data[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
        auto matrix = dense_CSR::Matrix_CSR(data);
        {
            // Timer t(out1);
            simple_iter::Jacobi_iter(matrix, vec, b, 0.0001, 3, 10);
            // auto res = simple_iter::Jacobi_iter(matrix, vec, b, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out2);
            simple_iter::Gauss_Seidel_iter(matrix, vec, b, 0.0001, 3, 10);
            // auto res = simple_iter::Gauss_Seidel_iter(matrix, vec, b, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out3);
            simple_iter::simple_iter(matrix, vec, b, tau_opt, 0.0001, 3, 10);
            // auto res = simple_iter::simple_iter(matrix, vec, b, tau_opt, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out4);
            simple_iter::simple_iter_boost(matrix, vec, b, lambda_min, lambda_max, 0.0001, 3, 10, 3);
            // auto res = simple_iter::simple_iter_boost(matrix, vec, b, lambda_min, lambda_max, 0.0001, 3, 100, 3).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }

}