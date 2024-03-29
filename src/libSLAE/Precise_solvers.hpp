#pragma once
#include "dense_CSR.hpp"

namespace Precise_solvers{
    using vector = std::vector<double>;

    vector run_through_3_diag_matrix(const vector& as, const vector& bs,
     const vector& cs, const vector& ds);
    
    vector solve_qr(const vector& b, const dense_CSR::Matrix& A);

    std::pair<dense_CSR::Matrix, dense_CSR::Matrix> QR(const dense_CSR::Matrix& A); 
    vector reflection(const vector& v, const vector& x);

}
