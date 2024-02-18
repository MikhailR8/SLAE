#include "dense_CSR.hpp"
#include <chrono>
#include <fstream>

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
    std::fstream out4("result4.txt", std::ios_base::app);
    for(auto k = 100u; k < 2000; k+=5){
        auto dense_array = std::vector<std::vector<double>>(k, std::vector<double>(k));
        auto csr_array = std::vector<std::vector<double>>(k, std::vector<double>(k));
        auto vec = std::vector<double>(k);
        for (auto i = 0u; i < k; i++){
            for(auto j = 0u; j < k; j++){
                std::srand(std::time(nullptr)); // use current time as seed for random generator
                int random_value = std::rand();
                dense_array[i][j] = static_cast<double>(random_value);
            }
        }
        for (auto i = 0u; i < k; i++){
            std::srand(std::time(nullptr)); // use current time as seed for random generator
            int random_value = std::rand();
            vec[i] = static_cast<double>(random_value);
        }
        for (auto i = 0u; i < k; i++){
            for(auto j = 0u; j < k; j++){
                if(i % 2 == 0 || j % 2 == 0){
                    csr_array[i][j] = 0;
                }
                else{
                    std::srand(std::time(nullptr)); // use current time as seed for random generator
                    int random_value = std::rand();
                    dense_array[i][j] = static_cast<double>(random_value);
                }
            }
        }
        auto matrix1 = cmake_example::Matrix(dense_array);
        auto matrix2 = cmake_example::Matrix_CSR(dense_array);
        auto matrix3 = cmake_example::Matrix(csr_array);
        auto matrix4 = cmake_example::Matrix_CSR(csr_array);
        {
            Timer t(out1);
            matrix1 * vec;
        }
        {
            Timer t(out2);
            matrix2 * vec;
        }
        {
            Timer t(out3);
            matrix3 * vec;
        }
        {
            Timer t(out4);
            matrix4 * vec;
        }
    }
    out1.close();
    out2.close();
    out3.close();
    out4.close();
}