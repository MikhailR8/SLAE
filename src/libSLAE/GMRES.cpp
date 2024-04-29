#include "GMRES.hpp"

namespace dense_CSR{
    GMRES_cache::GMRES_cache(const Matrix_CSR& matrix, unsigned max_iter):
     A{matrix}, vs{std::vector<vector>(matrix.get_size().first, vector(max_iter))}, 
     Givens_rotations{std::vector<vector>(max_iter, vector(2u))}, filled_rows{0u},
     Rs{std::vector<vector>(max_iter + 1u, vector(max_iter))} {}

    void GMRES_cache::add_v_and_rotation_and_R(const vector& inputv, double cos,
         double sin, const vector& inputR){
        vs.set_col(0u, inputv, filled_rows);
        Rs.set_col(0u, inputR, filled_rows);
        Givens_rotations.set_row(0u, std::vector({cos, sin}), filled_rows);
        filled_rows++;        
    }

    Matrix GMRES_cache::get_V() const{
        return vs;
    }

    vector GMRES_cache::get_v_num(unsigned num) const{
        return vs.get_col(0u, A.get_size().first, num);
    }

    vector GMRES_cache::get_rotation(unsigned num) const{
        return Givens_rotations.get_row(0u, 2u, num);
    }

    Matrix GMRES_cache::get_R() const{
        return Rs;
    }

    unsigned int GMRES_cache::get_filled() const{
        return filled_rows;
    }
}

namespace Iter_solvers{
    using dense_CSR::operator+; //Без этого из Iter_solvers не видно операторов из другого namespace'а
    using dense_CSR::operator*;
    using dense_CSR::operator-;

    vector GMRES(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned max_iteration, bool testmode){
        auto cache = Arnoldi_alg(A, x_0, b, max_iteration, target_discrepancy, testmode);
        unsigned filled = cache.get_filled();
        vector z(filled);
        auto answer = vector(x_0);
        test_pair out({vector(), std::vector<unsigned long long>()});

        if(testmode){
            {
                Timer T(&out.second);
                //z = Q^T * e^1 * ||r^0||, нас интересует только первый столбец Q^T = первая строка Q = 
                // = первый столбец произведения матриц поворотов R^l*...*R^1, прямым перемножением можно
                //получить, что этот столбец имеет вид:
                // (c_0, c_1*s_0, c_2*s_0*s_1, ..., c_n*s_0*...*s_(n-1), s_n*s_0*...*s_(n-1))^T, здесь 
                // c - cos и s - sin поворотов. Таким образом, последний элемент этого столбца гамма = 
                //произведению всех синусов поворотов, а все остальные элементы вычисляем ниже
                double sin_multiplication = 1.0;
                for(auto i = 0u; i < filled; i++){
                    auto temp = cache.get_rotation(i); //:( structured binding не работают с вектором
                    double cos = temp[0], sin = temp[1];
                    z[i] = cos * sin_multiplication;
                    sin_multiplication *= sin;
                }
                z = dense_CSR::get_length(A * x_0 - b) * z;
        
                auto R = cache.get_R();
                vector y(z);
                for (auto k = 0u; k < filled; k++){
                    auto i = filled - 1u - k; //новый счётчик
                    auto row = R.get_row(0, filled, i);
                    y[i] *= (1 / row[i]);
                    row = (1 / row[i]) * row;
                    R.set_row(0, row, i);
                    for (auto p = 0u; p < k; p++){
                        auto j = filled - 1u - p; //новый счётчик
                        y[i] -= y[j] * (R(i, j));
                    }
                }

                auto V = cache.get_V();
                for(auto i = 0u; i < answer.size(); i++){
                    answer[i] = answer[i] - (V.get_row(0, filled, i) * y);
                }
            }
            out.first.push_back(dense_CSR::get_length(A * answer - b));
        } else{
            //z = Q^T * e^1 * ||r^0||, нас интересует только первый столбец Q^T = первая строка Q = 
            // = первый столбец произведения матриц поворотов R^l*...*R^1, прямым перемножением можно
            //получить, что этот столбец имеет вид:
            // (c_0, c_1*s_0, c_2*s_0*s_1, ..., c_n*s_0*...*s_(n-1), s_n*s_0*...*s_(n-1))^T, здесь 
            // c - cos и s - sin поворотов. Таким образом, последний элемент этого столбца гамма = 
            //произведению всех синусов поворотов, а все остальные элементы вычисляем ниже
            double sin_multiplication = 1.0;
            for(auto i = 0u; i < filled; i++){
                auto temp = cache.get_rotation(i); //:( structured binding не работают с вектором
                double cos = temp[0], sin = temp[1];
                z[i] = cos * sin_multiplication;
                sin_multiplication *= sin;
            }
            z = dense_CSR::get_length(A * x_0 - b) * z;
    
            auto R = cache.get_R();
            vector y(z);
            for (auto k = 0u; k < filled; k++){
                auto i = filled - 1u - k; //новый счётчик
                auto row = R.get_row(0, filled, i);
                y[i] *= (1 / row[i]);
                row = (1 / row[i]) * row;
                R.set_row(0, row, i);
                for (auto p = 0u; p < k; p++){
                    auto j = filled - 1u - p; //новый счётчик
                    y[i] -= y[j] * (R(i, j));
                }
            }

            auto V = cache.get_V();
            for(auto i = 0u; i < answer.size(); i++){
                answer[i] = answer[i] - (V.get_row(0, filled, i) * y);
            }            
        }
        if(testmode) print_to_file(out, "GMRES");
        return answer;
    }

    dense_CSR::GMRES_cache Arnoldi_alg(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    unsigned max_iteration, double target_discrepancy, bool testmode){
        auto cache = dense_CSR::GMRES_cache(A, max_iteration);
        vector t(max_iteration);
        auto v = A * x_0 - b;
        double discrepancy_if_execute = dense_CSR::get_length(v); //gamma_0
        v = (1 / dense_CSR::get_length(v)) * v; //v_0
        test_pair out({vector(), std::vector<unsigned long long>()});

        for (auto j = 0u; j < max_iteration; j++){
            if(testmode){
                {
                    Timer T(&out.second);

                    vector h(max_iteration + 1u);
                    t = A * v;
                    for (auto k = 0u; k < j + 1u; k++){
                        if(k != j){
                            h[k] = cache.get_v_num(k) * t;
                            t = t - h[k] * cache.get_v_num(k);
                        }
                        else{
                            h[k] = v * t;
                            t = t - h[k] * v;                    
                        }
                    }
                    auto temp_h_j_plus_1 = dense_CSR::get_length(t);
                    h[j+1] = temp_h_j_plus_1;
                    //Вращения Гивенса: сначала применяем уже посчитанные, затем вычисляем новое
                    for(auto i = 0u; i < j; i++){
                        vector cos_sin = cache.get_rotation(i);
                        auto temp = h[i];
                        h[i] = temp * cos_sin[0] - h[i+1] * cos_sin[1];
                        h[i+1] = temp * cos_sin[1] + h[i+1] * cos_sin[0];
                    }
                    double new_cos = h[j] / std::sqrt(h[j] * h[j] + h[j+1] * h[j+1]);
                    double new_sin = (-1) * h[j+1] / std::sqrt(h[j] * h[j] + h[j+1] * h[j+1]);
                    double temp = h[j];
                    h[j] = temp * new_cos - h[j+1] * new_sin;
                    h[j+1] = temp * new_sin + h[j+1] * new_cos;
                    cache.add_v_and_rotation_and_R(v, new_cos, new_sin, h);

                    //Конец вращений, в кэш записаны v, новое вращение и приведённый к верхнетреугольному виду
                    //вектор h
                    v = (1 / temp_h_j_plus_1) * t;

                    discrepancy_if_execute *= new_sin; //См функцию GMRES, там пояснение, откуда это                    
                }
                out.first.push_back(discrepancy_if_execute);
                if(std::abs(discrepancy_if_execute) <= target_discrepancy){
                    print_to_file(out, "GMRES");
                    return cache;
                }
            } else{
                vector h(max_iteration + 1u);
                t = A * v;
                for (auto k = 0u; k < j + 1u; k++){
                    if(k != j){
                        h[k] = cache.get_v_num(k) * t;
                        t = t - h[k] * cache.get_v_num(k);
                    }
                    else{
                        h[k] = v * t;
                        t = t - h[k] * v;                    
                    }
                }
                auto temp_h_j_plus_1 = dense_CSR::get_length(t);
                h[j+1] = temp_h_j_plus_1;
                //Вращения Гивенса: сначала применяем уже посчитанные, затем вычисляем новое
                for(auto i = 0u; i < j; i++){
                    vector cos_sin = cache.get_rotation(i);
                    auto temp = h[i];
                    h[i] = temp * cos_sin[0] - h[i+1] * cos_sin[1];
                    h[i+1] = temp * cos_sin[1] + h[i+1] * cos_sin[0];
                }
                double new_cos = h[j] / std::sqrt(h[j] * h[j] + h[j+1] * h[j+1]);
                double new_sin = (-1) * h[j+1] / std::sqrt(h[j] * h[j] + h[j+1] * h[j+1]);
                double temp = h[j];
                h[j] = temp * new_cos - h[j+1] * new_sin;
                h[j+1] = temp * new_sin + h[j+1] * new_cos;
                cache.add_v_and_rotation_and_R(v, new_cos, new_sin, h);

                //Конец вращений, в кэш записаны v, новое вращение и приведённый к верхнетреугольному виду
                //вектор h
                v = (1 / temp_h_j_plus_1) * t;

                discrepancy_if_execute *= new_sin; //См функцию GMRES, там пояснение, откуда это
                if(std::abs(discrepancy_if_execute) <= target_discrepancy){
                    return cache;
                }
            }
        }
        if(testmode) print_to_file(out, "GMRES");
        return cache;
    }
}