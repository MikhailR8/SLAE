#include "Precise_solvers.hpp"

namespace Precise_solvers{

    vector run_through_3_diag_matrix(const vector& as, const vector& bs, const vector& cs, const vector& ds){
        vector::size_type n = as.size();
        vector ps = vector(n);
        vector qs = vector(n);
        vector xs = vector(n + 1);

        ps[0] = -cs[0] / bs[0];
        qs[0] = ds[0] / bs[0];

        for(auto i = 1u; i < n; i++){
            ps[i] = -cs[i] / (as[i - 1] * ps[i - 1] + bs[i]);
            qs[i] = (ds[i] - as[i - 1] * qs[i - 1]) / (as[i - 1] * ps[i - 1] + bs[i]);
        }

        xs[n] = (ds[n] - as[n - 1] * qs[n - 1]) / (ps[n - 1] * as[n - 1] + bs[n]);

        for(auto i = n - 1; i != 0; i--) {
            xs[i] = ps[i] * xs[i + 1] + qs[i];
        }
        xs[0] = ps[0] * xs[1] + qs[0];
        return xs;
    }

    vector solve_qr(const vector& b, const dense_CSR::Matrix& A){
        auto [Q, R] = QR(A);
        vector::size_type m = R.get_size().first;
        vector::size_type n = R.get_size().second;
        Q.transpose();
        auto answer = Q * b;
        for (auto k = 0u; k < n; k++){
            auto i = n - 1 - k; //новый счётчик
            auto row = R.get_row(0, n, i);
            answer[i] *= (1 / row[i]);
            row = dense_CSR::operator*((1 / row[i]), row);
            R.set_row(0, row, i);
            for (auto p = 0u; p < n - 1u - i; p++){
                auto j = n - 1 - p; //новый счётчик
                answer[i] -= answer[j] * (R(i, j));
            }
        }
        return answer;
    }
    
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

}