#include "Householder_alg.hpp"

namespace Householder{
    // using dense_CSR::operator*;
    // using dense_CSR::operator+;

    // void print_vector(const vector& input){
    //     for (auto i = 0u; i < input.size(); i++){
    //         std::cout << input[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    vector reflection(const vector& v, const vector& x){
        auto out = vector(x);
        auto koef = (-2.0) * dense_CSR::operator*(x, v) / dense_CSR::operator*(v, v);
        auto out1 = dense_CSR::operator*(koef, v);
        // print_vector(dense_CSR::operator+(out, out1));
        return dense_CSR::operator+(out, out1);
    }

    std::pair<dense_CSR::Matrix, dense_CSR::Matrix> QR(const dense_CSR::Matrix& A_input){
        auto A = dense_CSR::Matrix(A_input);
        vector::size_type m = A.get_size().first;
        vector::size_type n = A.get_size().second;
        auto to_P = std::vector<vector>(m, vector(m));
        for (auto i = 0u; i < m; i++){
            to_P[i][i] = 1.0;
        }
        auto P = dense_CSR::Matrix(to_P);
        auto v_arr = std::vector<vector>(n);
        for (auto i = 0u; i + 1 < n; i++){
            vector x = A.get_col(i, m, i);
            v_arr[i] = vector(x);
            if (x[0] >= 0.0){
                v_arr[i][0] += dense_CSR::get_length(x);
            }
            else{
                v_arr[i][0] -= dense_CSR::get_length(x);
            }

            // std::cout << i << std::endl;
            // print_vector(v_arr[i]);
            // dense_CSR::print_matrix(A);
            // std::cout << " " << std::endl;
            // dense_CSR::print_matrix(P);

            for (auto j = i; j < n; j++){
                vector col_j = A.get_col(i, m, j);
                A.set_col(i, reflection(v_arr[i], col_j), j);
            }
            for (auto j = 0u; j < m; j++){
                vector row_j = P.get_row(i, m, j);
                P.set_row(i, reflection(v_arr[i], row_j), j);
            }

        }
        return std::pair<dense_CSR::Matrix, dense_CSR::Matrix>(P, A);
    }

    vector solve_qr(const vector& b, const dense_CSR::Matrix& A){
        auto pair = QR(A);
        auto Q = pair.first;
        auto R = pair.second;
        vector::size_type m = R.get_size().first;
        vector::size_type n = R.get_size().second;
        Q.transpose();
        auto answer = Q * b;
        for (int i = static_cast<int>(b.size()) - 1; i >= 0; i--){
            auto row = R.get_row(0, n, i);
            answer[i] *= (1 / row[i]);
            row = dense_CSR::operator*((1 / row[i]), row);
            R.set_row(0, row, i);
            for (int j = static_cast<int>(b.size()) - 1; j > i; j--){
                answer[i] -= answer[j] * (R(i, j));
            }
        }
        return answer;
    }

}