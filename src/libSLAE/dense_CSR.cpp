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

    vector Matrix_CSR::Jacobi_iteration(const vector& v, const vector& b) const{
        auto res = vector(b);
        for (auto i = 0u; i < v.size(); ++i) {
            auto divider = 0.0;
            for (auto k = rows[i]; k < rows[i + 1]; ++k) {
                if(i != cols[k]){
                    res[i] -= values[k] * v[cols[k]];
                }
                else{
                    divider = values[k];
                }
            }
            res[i] /= divider;
        }
        return res;
    }

    vector Matrix_CSR::Gauss_Seidel_iteration(const vector& v, const vector& b) const{
        auto res = vector(v);
    //x_{0}^{i+1} = 1/D_0 * (b_0 - A_{01}x_{1}^{i} - A_{02}x_{2}^{i} - ... - A_{0, n-1}x_{n-1}^{i})
    //x_{1}^{i+1} = 1/D_1 * (b_1 - A_{12}x_{2}^{i} - ... - A_{1, n-1}x_{n-1}^{i} - A_{10}x_{0}^{i+1})
    //...
    //Последнее слагаемое в выражении x_{1}^{i+1} содержит уже вычисленное x_{0}^{i+1}, а при
    //вычислении x{j}^{i+1} элемент x{j}^{i} не используется, поэтому в начале приравниваем
    //res[i] = b[i] (это безопасно), а проходя по CSR матрице запоминаем диагональный элемент, чтобы
    //сэкономить время
        for(auto i = 0u; i < n; i++){
            res[i] = b[i];
            auto divider = 0.0;
            for (auto k = rows[i]; k < rows[i + 1]; ++k) {
                if(cols[k] != i){
                    res[i] -= values[k] * res[cols[k]];
                }
                else{
                    divider = values[k];
                }
            }
            res[i] /= divider;
        }
        return res;
    }

    vector Matrix_CSR::Symmetric_Gauss_Seidel_iteration(const vector& v, const vector& b) const {
        auto res = vector(v);
    //x_{0}^{i+1} = 1/D_0 * (b_0 - A_{01}x_{1}^{i} - A_{02}x_{2}^{i} - ... - A_{0, n-1}x_{n-1}^{i})
    //x_{1}^{i+1} = 1/D_1 * (b_1 - A_{12}x_{2}^{i} - ... - A_{1, n-1}x_{n-1}^{i} - A_{10}x_{0}^{i+1})
    //...
    //Последнее слагаемое в выражении x_{1}^{i+1} содержит уже вычисленное x_{0}^{i+1}, а при
    //вычислении x{j}^{i+1} элемент x{j}^{i} не используется, поэтому в начале приравниваем
    //res[i] = b[i] (это безопасно), а проходя по CSR матрице запоминаем диагональный элемент, чтобы
    //сэкономить время
        for(auto i = 0u; i < n; i++){
            res[i] = b[i];
            auto divider = 0.0;
            for (auto k = rows[i]; k < rows[i + 1]; ++k) {
                if(cols[k] != i){
                    res[i] -= values[k] * res[cols[k]];
                }
                else{
                    divider = values[k];
                }
            }
            res[i] /= divider;
        }

        for(auto i = 0u; i < n; i++){
            unsigned j = n - 1 - i; //новый счётчик
            res[j] = b[j];
            auto divider = 0.0;
            for (auto k = rows[j]; k < rows[j + 1]; ++k) {
                if(cols[k] != j){
                    res[j] -= values[k] * res[cols[k]];
                }
                else{
                    divider = values[k];
                }
            }
            res[j] /= divider;
        }
        return res;        
    }

    double Matrix_CSR::power_iteration(unsigned max_iteration, double precision) const {
        vector r(n);
        for(auto i = 0u; i < n; i++){
            r[i] = 1.0; //r_0
        }
        vector Ar = this->operator*(r); //Ar_0
        r = (1 / get_length(Ar)) * Ar; //r_1
        Ar = this->operator*(r); //Ar_1
        double max_eigenvalue = (r * Ar) / (r * r); // mu_1
        for(auto i = 1u; i < max_iteration; i++){
            r = (1 / get_length(Ar)) * Ar;
            Ar = this->operator*(r);
            double temp_max_eigenvalue = (r * Ar) / (r * r);
            if(std::abs(temp_max_eigenvalue - max_eigenvalue) <= precision){
                return temp_max_eigenvalue;
            }
            else{
                max_eigenvalue = temp_max_eigenvalue;
            }
        }
        return max_eigenvalue;
    }


    GMRES_cache::GMRES_cache(const Matrix_CSR& matrix, unsigned max_iter):
     A{matrix}, vs{std::vector<vector>(matrix.get_size().first, vector(max_iter))}, 
     Givens_rotations{std::vector<vector>(max_iter, vector(2u))}, filled_rows{0u},
     Rs{std::vector<vector>(matrix.get_size().first, vector(max_iter))} {}

    void GMRES_cache::add_v_and_rotation_and_R(const vector& inputv, double cos,
         double sin, const vector& inputR){
        vs.set_col(0u, inputv, filled_rows);
        Rs.set_col(0u, inputR, filled_rows);
        Givens_rotations.set_row(0u, std::vector({cos, sin}), filled_rows);
        filled_rows++;        
    }

    Matrix GMRES_cache::get_V(){
        return vs;
    }

    vector GMRES_cache::get_rotation(unsigned num){
        return Givens_rotations.get_row(0u, 2u, num);
    }

    Matrix GMRES_cache::get_R(){
        return Rs;
    }


    vector operator+(const vector& lhs, const vector& rhs){
        auto res = vector(lhs.size());
        for (auto i = 0u; i < lhs.size(); i++){
            res[i] = lhs[i] + rhs[i];
        }
        return res;
    }

    vector operator-(const vector& lhs, const vector& rhs){
        auto res = vector(lhs.size());
        for (auto i = 0u; i < lhs.size(); i++){
            res[i] = lhs[i] - rhs[i];
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

    void print_vector(const vector& input){
        for (auto i = 0u; i < input.size(); i++){
            std::cout << input[i] << " ";
        }
        std::cout << std::endl;
    }

}