#include "dense_CSR.hpp"
#include "simple_iter.hpp"
#include <chrono>
#include <fstream>
#include <random>

using steady_clock = std::chrono::steady_clock;
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

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist0(0,10);
    std::uniform_int_distribution<std::mt19937::result_type> dist1(20000,30000);
    std::uniform_int_distribution<std::mt19937::result_type> dist2(1,10);

    //Прогоняем матрицы разных размеров
    for(auto k = 1000u; k < 2000u; k+=20){
        auto dense_array = std::vector<std::vector<double>>(k, std::vector<double>(k));
        auto vec = std::vector<double>(k);
        auto b = std::vector<double>(k);
        //Заполняем матрицу плотно
        for (auto i = 0u; i < k; i++){
            for(auto j = 0u; j < k; j++){
                if(i != j){
                    dense_array[i][j] = static_cast<double>(dist0(rng));
                }
                else{
                    dense_array[i][j] = static_cast<double>(dist1(rng));
                }
            }
        }
        //Заполняем вектор
        for (auto i = 0u; i < k; i++){
            vec[i] = static_cast<double>(dist2(rng));
        }

        for (auto i = 0u; i < k; i++){
            b[i] = static_cast<double>(dist2(rng));
        }
        
        auto matrix = dense_CSR::Matrix_CSR(dense_array);
        {
            Timer t(out1);
            simple_iter::Jacobi_iter(matrix, vec, b, 0.01, 3, 100);
        }
        {
            Timer t(out2);
            simple_iter::Gauss_Seidel_iter(matrix, vec, b, 0.01, 3, 100);
        }
        {
            Timer t(out3);
            simple_iter::simple_iter(matrix, vec, b, 0.12, 0.01, 3, 100);
        }
    }
    out1.close();
    out2.close();
    out3.close();
}