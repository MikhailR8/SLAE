#include "dense_CSR.hpp"
#include "simple_iter.hpp"
#include <chrono>
#include <fstream>
#include <random>
#include <iostream>

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

int main(){
    std::fstream out1("result1.txt", std::ios_base::app);
    std::fstream out2("result2.txt", std::ios_base::app);
    std::fstream out3("result3.txt", std::ios_base::app);
    std::fstream out4("result4.txt", std::ios_base::app);

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
    
    vector res1(100);
    vector res2(100);
    vector res3(100);
    vector res4(100);
        auto matrix = dense_CSR::Matrix_CSR(data);
        {
            // Timer t(out1);
            res1 = simple_iter::Jacobi_iter(matrix, vec, b, 0.0001, 3, 100);
            // auto res = simple_iter::Jacobi_iter(matrix, vec, b, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out2);
            res2 = simple_iter::Gauss_Seidel_iter(matrix, vec, b, 0.0001, 3, 100);
            // auto res = simple_iter::Gauss_Seidel_iter(matrix, vec, b, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out3);
            res3 = simple_iter::simple_iter(matrix, vec, b, tau_opt, 0.0001, 3, 100);
            // auto res = simple_iter::simple_iter(matrix, vec, b, tau_opt, 0.0001, 3, 100).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
        {
            // Timer t(out4);
            res4 = simple_iter::simple_iter_boost(matrix, vec, b, lambda_min, lambda_max, 0.0001, 3, 100, 3);
            // auto res = simple_iter::simple_iter_boost(matrix, vec, b, lambda_min, lambda_max, 0.0001, 3, 100, 3).first;
            // for(auto i = 0u; i < n; i++){
            //     std::cout << res[i] << std::endl;
            // }
            // std::cout << std::endl;
        }
    for(auto i = 0u; i < 100u; i++){
            out1 << res1[i] << std::endl;
            out2 << res2[i] << std::endl;
            out3 << res3[i] << std::endl;
            out4 << res4[i] << std::endl;
        }
    out1.close();
    out2.close();
    out3.close();
    out4.close();
}