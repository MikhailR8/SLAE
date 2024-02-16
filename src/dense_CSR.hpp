#pragma once
#include <vector>
#include <iostream>

namespace cmake_example{
    using vector = std::vector<double>;
    class Matrix{
    public:
        Matrix(std::vector<std::vector<double>>& input);

        double operator()(int i, int j) const;

        std::pair<unsigned, unsigned> get_size() const;

        vector operator*(const vector &v) const; 

    private:
        vector data;
        unsigned int m;
        unsigned int n;
    };

    class Matrix_CSR{
    public:
        Matrix_CSR(std::vector<std::vector<double>>& input);

        double operator()(int i, int j) const;

        std::pair<unsigned, unsigned> get_size() const;

        vector operator*(const vector& v) const;
    private:
        vector values;
        vector cols;
        vector rows;
        unsigned int m;
        unsigned int n;
    };

    template <typename T>
    void print_matrix(T matrix){
        std::pair<unsigned, unsigned> size_matrix = matrix.get_size();
        for(auto i = 0u; i < size_matrix.first; i++){
            for(auto j = 0u; j < size_matrix.second; j++){
                std::cout << matrix(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    vector operator+(const vector& lhs, const vector& rhs);

    double operator*(const vector& lhs, const vector& rhs);
}