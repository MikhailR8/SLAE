#include <iostream>
#include <vector>

using vector = std::vector<double>;

vector run_through_3_diag_matrix(vector as, vector bs, vector cs, vector ds){
    auto n = as.size();
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
