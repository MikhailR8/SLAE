#pragma once
#include <fstream>
#include <chrono>
#include <vector>
#include <iostream>

using steady_clock = std::chrono::steady_clock;
using vector = std::vector<double>;
using test_pair = std::pair<vector, std::vector<unsigned long long>>;
class Timer {
public:
    Timer(std::vector<unsigned long long>* to_write) : start(steady_clock::now()), points(to_write) {}
    ~Timer() {
        //int_64
        auto point = std::chrono::duration_cast<std::chrono::microseconds>(steady_clock::now()
            - start).count();
        points->push_back(static_cast<unsigned long long>(point));
    }
private:
    steady_clock::time_point start;
    //Да, указатели не лучшая идея, но это будет слишком долго, если мы будем копировать вектор как поле
    //класса каждый раз измеряя время контрольной точки, а иначе его не сохранить как поле класса
    std::vector<unsigned long long>* points;
};

void print_to_file(const test_pair& to_write, std::string filename);