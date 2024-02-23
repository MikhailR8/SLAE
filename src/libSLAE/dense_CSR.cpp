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

    Matrix::Matrix(const Matrix& input): m{input.get_size().first}, n{input.get_size().second},
    data{input.get_raw_data()} {}

    double Matrix::operator()(unsigned i, unsigned j) const{
        return data[n * i + j];
    }

    void Matrix::set_ij(unsigned i, unsigned j, double value){
        data[n * i + j] = value;
    }

    void Matrix::set_col(unsigned start_point, vector col, unsigned num_col){
        for(auto i = 0u; i < col.size(); i++){
            data[n * (i + start_point) + num_col] = col[i];
        }
    }

    void Matrix::set_row(unsigned start_point, vector row, unsigned num_row){
        for(auto i = 0; i < row.size(); i++){
            data[n * num_row + (i + start_point)] = row[i];
        }
    }

    vector Matrix::get_col(unsigned start_point, unsigned stop_point, unsigned num_col) const{
        int size = static_cast<int>(stop_point) - static_cast<int>(start_point);
        auto out = vector(size);
        for(auto i = 0; i < size; i++){
            out[i] = data[n * (i + static_cast<int>(start_point)) + num_col];
        }
        return out;
    }

    vector Matrix::get_row(unsigned start_point, unsigned stop_point, unsigned num_row) const{
        int size = static_cast<int>(stop_point) - static_cast<int>(start_point);
        auto out = vector(size);
        for(auto i = 0; i < size; i++){
            out[i] = data[n * num_row + (i + static_cast<int>(start_point))];
        }
        return out;
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

    vector Matrix::get_raw_data() const{
        return vector(data);
    }  

    void Matrix::transpose(){
        for (auto i = 0u; i < m; i++){
            for(auto j = i + 1; j < n; j++){
                double temp = data[n * i + j];
                data[n * i + j] = data[n * j + i];
                data[n * j + i] = temp;
            }
        }
        auto temp = n;
        n = m;
        m = temp;
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
            res += lhs[i] * rhs[i];
        }
        return res;       
    }

    vector operator*(double rhs, const vector& lhs){
        vector out = vector(lhs.size());
        for (auto i = 0u; i < lhs.size(); i++){
            out[i] = rhs * lhs[i];
        }
        return out;
    }

    double get_length(const vector& v){
        double out = 0;
        for (auto i = 0u; i < v.size(); i++){
            out += (v[i] * v[i]);
        }
        return std::sqrt(out);
    }

}