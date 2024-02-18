#include "dense_CSR.hpp"

namespace dense_CSR{
    Matrix::Matrix(const std::vector<vector>& input): m{input.size()}, n{input[0].size()},
    data{vector(input.size() * input[0].size())}{
        for(auto i = 0u; i < m; i++){
            for(auto j = 0u; j < n; j++){
                data[n * i + j] = input[i][j];
            }
        }
    }

    double Matrix::operator()(int i, int j) const{
        return data[n * i + j];
    }

    std::pair<vector::size_type, vector::size_type> Matrix::get_size() const{
        return std::pair<vector::size_type, vector::size_type>({m, n});
    }

    vector Matrix::operator*(const vector &v) const{
        vector res(v.size());
        for (auto i = 0u; i < m; i++) {
            for(auto j = 0u; j < n; j++){
                res[i] += data[n * i + j] * v[j];
            }
        }
        return res;
    }        

    Matrix_CSR::Matrix_CSR(const std::vector<vector>& input): m{input.size()}, n{input[0].size()},
    values{vector()}, cols{std::vector<unsigned>()}, rows{std::vector<unsigned>({0u})}{
        auto counter_elems = 0u;
        for(auto i = 0u; i < m; i++){
            for(auto j = 0u; j < n; j++){
                if(input[i][j] != 0){
                    values.push_back(input[i][j]);
                    cols.push_back(j);
                    counter_elems++;
                }
            }
            rows.push_back(counter_elems);
        }
    }

    double Matrix_CSR::operator()(int i, int j) const{
        for (unsigned k = rows[i]; k < rows[i + 1]; ++k) {
            if (cols[k] == j) {
                return values[k];
            }
        }
        return 0;
    }

    std::pair<vector::size_type, vector::size_type> Matrix_CSR::get_size() const{
        return std::pair<vector::size_type, vector::size_type>({m, n});
    }

    vector Matrix_CSR::operator*(const vector& v) const {
        vector res(v.size());
        for (int i = 0; i < v.size(); ++i) {
            for (int k = this->rows[i]; k < this->rows[i + 1]; ++k) {
                res[i] += this->values[k] * v[this->cols[k]];
            }
        }
        return res;
    }


    vector operator+(const vector& lhs, const vector& rhs){
        auto res = vector(lhs.size());
        for (auto i = 0u; i < lhs.size(); i++){
            res[i] = lhs[i] + rhs[i];
        }
        return res;
    }

    double operator*(const vector& lhs, const vector& rhs){
        double res = 0;
        for (auto i = 0u; i < lhs.size(); i++){
            res += lhs[i] + rhs[i];
        }
        return res;       
    }

}