#include "dense_CSR.hpp"

namespace cmake_example{
    Matrix::Matrix(std::vector<std::vector<double>>& input){
        unsigned int m_input = input.size();
        unsigned int n_input = input[0].size();
        m = m_input;
        n = n_input;
        data = vector(m * n);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                data[n * i + j] = input[i][j];
            }
        }
    }

    double Matrix::operator()(int i, int j) const{
        return data[n * i + j];
    }

    std::pair<unsigned, unsigned> Matrix::get_size() const{
        return std::pair<unsigned, unsigned>({m, n});
    }

    vector Matrix::operator*(const vector &v) const{
        std::vector<double> res(v.size());
        for (auto i = 0u; i < m; i++) {
            for(auto j = 0u; j < n; j++){
                res[i] += data[n * i + j] * v[j];
            }
        }
        return res;
    }        

    Matrix_CSR::Matrix_CSR(std::vector<std::vector<double>>& input){
        unsigned int m_input = input.size();
        unsigned int n_input = input[0].size();
        m = m_input;
        n = n_input;
        values = vector();
        cols = vector();
        rows = vector({0});
        int counter_elems = 0;
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
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
        for (int k = rows[i]; k < rows[i + 1]; ++k) {
            if (cols[k] == j) {
                return values[k];
            }
        }
        return 0;
    }

    std::pair<unsigned, unsigned> Matrix_CSR::get_size() const{
        return std::pair<unsigned, unsigned>({m, n});
    }

    vector Matrix_CSR::operator*(const vector& v) const {
        std::vector<double> res(v.size());
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