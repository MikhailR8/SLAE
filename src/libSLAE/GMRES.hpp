#pragma once
#include "dense_CSR.hpp"
#include "Test_speed_tools.hpp"

namespace dense_CSR{
    class GMRES_cache{
    public:
        GMRES_cache(const Matrix_CSR& matrix, unsigned max_iter);

        void add_v_and_rotation_and_R(const vector& inputv, double cos,
         double sin, const vector& inputR);

        Matrix get_V() const;

        vector get_rotation(unsigned num) const;

        Matrix get_R() const;

        vector get_v_num(unsigned num) const;

        unsigned int get_filled() const;

        double get_dis() const;

        void set_dis(double dis);
    private:
        Matrix_CSR A; //Матрица СЛАУ n x n
        Matrix vs; //Матрица из векторов v, хранящихся в столбцы, имеет размер n x k, k - число итераций
        Matrix Givens_rotations; //Матрица k x 2, строка этой матрицы - (cos \theta_i, sin \theta_i)
        Matrix Rs; //Верхнетреугольная матрица, преобразованная из матрицы Хессинберга, (k+1) x k
        unsigned int filled_rows; //Количество заполненных строк в матрицe Givens_rotation/столбцов vs
        double discrepancy_if_exe;
    };
}

namespace Iter_solvers{
    vector GMRES(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    double target_discrepancy, unsigned max_iteration, bool testmode = false);

    dense_CSR::GMRES_cache Arnoldi_alg(const dense_CSR::Matrix_CSR& A, const vector& x_0, const vector& b,
    unsigned max_iteration, double target_discrepancy, bool testmode = false);
}