#pragma once
#include <vector>
#include "dense_CSR.hpp"

namespace Householder{
    using vector = std::vector<double>;
    std::pair<dense_CSR::Matrix, dense_CSR::Matrix> QR(const dense_CSR::Matrix& A); 

    vector reflection(const vector& v, const vector& x);
    vector solve_qr(const vector& b, const dense_CSR::Matrix& A);
}
