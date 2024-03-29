#pragma once
#include <vector>
#include <iostream>
#include <cmath>

namespace dense_CSR{
    using vector = std::vector<double>;
    class Matrix{
    public:
        Matrix(const std::vector<vector>& input);

        Matrix(const Matrix& input);

        double operator()(unsigned i, unsigned j) const;

        void set_ij(unsigned i, unsigned j, double value);

        void set_col(unsigned start_point, vector col, unsigned num_col);

        void set_row(unsigned start_point, vector row, unsigned num_row);

        vector get_col(unsigned start_point, unsigned stop_point, unsigned num_col) const;

        vector get_row(unsigned start_point, unsigned stop_point, unsigned num_row) const;

        std::pair<vector::size_type, vector::size_type> get_size() const;

        vector operator*(const vector &v) const; 

        vector get_raw_data() const;

        void transpose();

    private:
        vector data;
        vector::size_type m;
        vector::size_type n;
    };

    class Matrix_CSR{
    public:
        Matrix_CSR(const std::vector<vector>& input);

        double operator()(int i, int j) const;

        std::pair<vector::size_type, vector::size_type> get_size() const;

        vector operator*(const vector& v) const;

        vector Jacobi_iteration(const vector& v, const vector& b) const;

        vector Gauss_Seidel_iteration(const vector& v, const vector& b) const;

        vector Symmetric_Gauss_Seidel_iteration(const vector& v, const vector& b) const;
    private:
        vector values;
        std::vector<unsigned> cols;
        std::vector<unsigned> rows;
        vector::size_type m;
        vector::size_type n;
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

    vector operator-(const vector& lhs, const vector& rhs);

    double operator*(const vector& lhs, const vector& rhs);

    vector operator*(double rhs, const vector& lhs);

    double get_length(const vector& v);

    void print_vector(const vector& input);
}