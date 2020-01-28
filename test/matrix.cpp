#include "lemon/matrix.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

static bool roughly(double a, double b, double tol = 1e-4) {
    return std::fabs(a - b) < tol;
}

TEST_CASE("Cardano") {
    // Simple case where the answer is obvious
    auto simple = lemon::cardano(27, 0, 0, 0);

    CHECK(roughly(std::get<0>(simple).real(), 0.0));
    CHECK(roughly(std::get<1>(simple).real(), 0.0));
    CHECK(roughly(std::get<2>(simple).real(), 0.0));

    // Real roots
    auto real = lemon::cardano(1, 0, -15, -4);
    CHECK(roughly(std::get<0>(real).real(), 4.0));
    CHECK(roughly(std::get<1>(real).real(), -2.0 - std::sqrt(3.0)));
    CHECK(roughly(std::get<2>(real).real(), std::sqrt(3.0) - 2.0));
    CHECK(roughly(std::get<0>(real).imag(), 0.0));
    CHECK(roughly(std::get<1>(real).imag(), 0.0));
    CHECK(roughly(std::get<2>(real).imag(), 0.0));

    // Two roots are the same
    auto same = lemon::cardano(1, 7, 0, 0);
    CHECK(roughly(std::get<0>(same).real(), 0.0));
    CHECK(roughly(std::get<1>(same).real(), -7.0));
    CHECK(roughly(std::get<2>(same).real(), 0.0));
    CHECK(roughly(std::get<0>(same).imag(), 0.0));
    CHECK(roughly(std::get<1>(same).imag(), 0.0));
    CHECK(roughly(std::get<2>(same).imag(), 0.0));

    // Similar, to above. Used to ensure that the eps value is appropriate
    auto same2 = lemon::cardano(1, -7, 0, 0);
    CHECK(roughly(std::get<0>(same2).real(), 7.0));
    CHECK(roughly(std::get<1>(same2).real(), 0.0));
    CHECK(roughly(std::get<2>(same2).real(), 0.0));
    CHECK(roughly(std::get<0>(same2).imag(), 0.0));
    CHECK(roughly(std::get<1>(same2).imag(), 0.0));
    CHECK(roughly(std::get<2>(same2).imag(), 0.0));

    auto complex_roots = lemon::cardano(1, 2, 3, 4);

    CHECK(roughly(std::get<0>(complex_roots).real(), -1.6506));
    CHECK(roughly(std::get<1>(complex_roots).real(), -0.1747));
    CHECK(roughly(std::get<2>(complex_roots).real(), -0.1747));
    CHECK(roughly(std::get<0>(complex_roots).imag(), 0.0));
    CHECK(roughly(std::get<1>(complex_roots).imag(), 1.5469));
    CHECK(roughly(std::get<2>(complex_roots).imag(), -1.5469));
}

TEST_CASE("Eigenvalues and Eigenvectors") {
    auto matrix_1 = lemon::Matrix3D(1, -3, 3, 3, -5, 3, 6, -6, 4);
    auto eigenvalues_1 = lemon::eigenvalues(matrix_1);
    CHECK(roughly(std::get<0>(eigenvalues_1).real(), -2));
    CHECK(roughly(std::get<1>(eigenvalues_1).real(), -2));
    CHECK(roughly(std::get<2>(eigenvalues_1).real(), 4));
    CHECK(roughly(std::get<0>(eigenvalues_1).imag(), 0.0));
    CHECK(roughly(std::get<1>(eigenvalues_1).imag(), 0.0));
    CHECK(roughly(std::get<2>(eigenvalues_1).imag(), 0.0));

    // just to be safe, another easy one
    auto matrix_2 = lemon::Matrix3D(2, 0, 0, 0, 3, 4, 0, 4, 9);
    auto eigenvalues_2 = lemon::eigenvalues(matrix_2);
    CHECK(roughly(std::get<0>(eigenvalues_2).real(), 1));
    CHECK(roughly(std::get<1>(eigenvalues_2).real(), 2));
    CHECK(roughly(std::get<2>(eigenvalues_2).real(), 11));
    CHECK(roughly(std::get<0>(eigenvalues_2).imag(), 0.0));
    CHECK(roughly(std::get<1>(eigenvalues_2).imag(), 0.0));
    CHECK(roughly(std::get<2>(eigenvalues_2).imag(), 0.0));

    auto matrix_3 = lemon::Matrix3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
    auto eigenvalues_3 = lemon::eigenvalues(matrix_3);
    CHECK(roughly(std::get<0>(eigenvalues_3).real(), -0.5));
    CHECK(roughly(std::get<1>(eigenvalues_3).real(), -0.5));
    CHECK(roughly(std::get<2>(eigenvalues_3).real(), 1.0));
    CHECK(roughly(std::get<0>(eigenvalues_3).imag(), std::sqrt(3) / 2));
    CHECK(roughly(std::get<1>(eigenvalues_3).imag(), -std::sqrt(3) / 2));
    CHECK(roughly(std::get<2>(eigenvalues_3).imag(), 0.0));

    // Eigenvectors

    // Matrix 2 is actually the easiest case, so we start here

    auto eigenvectors_2 = lemon::eigenvectors(matrix_2, eigenvalues_2);
    CHECK(roughly(eigenvectors_2[0][0], 0.0 / std::sqrt(5.0)));
    CHECK(roughly(eigenvectors_2[0][1], 2.0 / std::sqrt(5.0)));
    CHECK(roughly(eigenvectors_2[0][2], -1.0 / std::sqrt(5.0)));

    CHECK(roughly(eigenvectors_2[1][0], -1.0 / std::sqrt(1.0)));
    CHECK(roughly(eigenvectors_2[1][1], 0.0 / std::sqrt(1.0)));
    CHECK(roughly(eigenvectors_2[1][2], 0.0 / std::sqrt(1.0)));

    CHECK(roughly(eigenvectors_2[2][0], 0.0 / std::sqrt(5.0)));
    CHECK(roughly(eigenvectors_2[2][1], 1.0 / std::sqrt(5.0)));
    CHECK(roughly(eigenvectors_2[2][2], 2.0 / std::sqrt(5.0)));

    auto wikipedia_example = lemon::Matrix3D(3, 2, 6, 2, 2, 5, -2, -1, -4);
    auto eigenvalues_wkpd = lemon::eigenvalues(wikipedia_example);
    auto eigenvectors_wkpd = lemon::eigenvectors(wikipedia_example, eigenvalues_wkpd);

    CHECK(roughly(eigenvectors_wkpd[0][0], -2.0 / std::sqrt(6.0)));
    CHECK(roughly(eigenvectors_wkpd[0][1], -1.0 / std::sqrt(6.0)));
    CHECK(roughly(eigenvectors_wkpd[0][2],  1.0 / std::sqrt(6.0)));

    CHECK(roughly(eigenvectors_wkpd[1][0], -1.0 / std::sqrt(3.0)));
    CHECK(roughly(eigenvectors_wkpd[1][1], -1.0 / std::sqrt(3.0)));
    CHECK(roughly(eigenvectors_wkpd[1][2], 1.0 / std::sqrt(3.0)));

    CHECK(roughly(eigenvectors_wkpd[2][0], 2.0 / std::sqrt(6.0)));
    CHECK(roughly(eigenvectors_wkpd[2][1], 1.0 / std::sqrt(6.0)));
    CHECK(roughly(eigenvectors_wkpd[2][2], -1.0 / std::sqrt(6.0)));
}
