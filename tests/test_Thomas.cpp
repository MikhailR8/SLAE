#include <gtest/gtest.h>
#include <libSLAE/Thomas_alg.hpp>

using vector = std::vector<double>;
TEST(Matrix_test, first_matrix) {
    vector as = {23, 31, 15};
    vector bs = {32, 34, 76, 54};
    vector cs = {12, 54, 13};
    vector ds = {122, 324, 127, 64};
    vector result = Thomas::run_through_3_diag_matrix(as, bs, cs, ds);
    ASSERT_NEAR(-21.6998, result[0], 0.0001);
    ASSERT_NEAR(68.0329, result[1], 0.0001);
    ASSERT_NEAR(-27.5930, result[2], 0.0001);
    ASSERT_NEAR(8.8499, result[3], 0.0001);
}

TEST(Matrix_test, second_matrix) {
    vector as = {2, 3};
    vector bs = {10, 20, 10};
    vector cs = {2, 3};
    vector ds = {5, 4, 7};
    vector result = Thomas::run_through_3_diag_matrix(as, bs, cs, ds);
    ASSERT_NEAR(0.4904, result[0], 0.0001);
    ASSERT_NEAR(0.0481, result[1], 0.0001);
    ASSERT_NEAR(0.6856, result[2], 0.0001);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}