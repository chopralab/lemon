#include "lemon/matrix.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using Catch::Detail::Approx;

TEST_CASE("Cardano") {
    // Simple case where the answer is obvious
    auto simple = lemon::cardano(27, 0, 0, 0);

    CHECK(std::get<0>(simple).real() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(simple).real() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(simple).real() == Approx(0.0).epsilon(1e-4));

    // Real roots
    auto real = lemon::cardano(1, 0, -15, -4);
    CHECK(std::get<0>(real).real() == Approx(4.0).epsilon(1e-4));
    CHECK(std::get<1>(real).real() == Approx(-2.0 - std::sqrt(3.0)).epsilon(1e-4));
    CHECK(std::get<2>(real).real() == Approx(std::sqrt(3.0) - 2.0).epsilon(1e-4));
    CHECK(std::get<0>(real).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(real).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(real).imag() == Approx(0.0).epsilon(1e-4));

    // Two roots are the same
    auto same = lemon::cardano(1, 7, 0, 0);
    CHECK(std::get<0>(same).real() == Approx(0.0).epsilon(1e-4).epsilon(1e-4));
    CHECK(std::get<1>(same).real() == Approx(-7.0).epsilon(1e-4));
    CHECK(std::get<2>(same).real() == Approx(0.0).epsilon(1e-4).epsilon(1e-4));
    CHECK(std::get<0>(same).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(same).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(same).imag() == Approx(0.0).epsilon(1e-4));

    // Similar, to above. Used to ensure that the eps value is appropriate
    auto same2 = lemon::cardano(1, -7, 0, 0);
    CHECK(std::get<0>(same2).real() == Approx(7.0).epsilon(1e-4));
    CHECK(std::get<1>(same2).real() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(same2).real() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<0>(same2).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(same2).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(same2).imag() == Approx(0.0).epsilon(1e-4));

    auto complex_roots = lemon::cardano(1, 2, 3, 4);

    CHECK(std::get<0>(complex_roots).real() == Approx(-1.6506).epsilon(1e-4));
    CHECK(std::get<1>(complex_roots).real() == Approx(-0.1747).epsilon(1e-4));
    CHECK(std::get<2>(complex_roots).real() == Approx(-0.1747).epsilon(1e-4));
    CHECK(std::get<0>(complex_roots).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(complex_roots).imag() == Approx(1.5469).epsilon(1e-4));
    CHECK(std::get<2>(complex_roots).imag() == Approx(-1.5469).epsilon(1e-4));
}

TEST_CASE("Eigenvalues and Eigenvectors") {
    auto matrix_1 = lemon::Matrix3D(1, -3, 3, 3, -5, 3, 6, -6, 4);
    auto eigenvalues_1 = lemon::eigenvalues(matrix_1);
    CHECK(std::get<0>(eigenvalues_1).real() == Approx(-2).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_1).real() == Approx(-2).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_1).real() == Approx(4).epsilon(1e-4));
    CHECK(std::get<0>(eigenvalues_1).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_1).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_1).imag() == Approx(0.0).epsilon(1e-4));

    // just to be safe, another easy one
    auto matrix_2 = lemon::Matrix3D(2, 0, 0, 0, 3, 4, 0, 4, 9);
    auto eigenvalues_2 = lemon::eigenvalues(matrix_2);
    CHECK(std::get<0>(eigenvalues_2).real() == Approx(1).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_2).real() == Approx(2).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_2).real() == Approx(11).epsilon(1e-4));
    CHECK(std::get<0>(eigenvalues_2).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_2).imag() == Approx(0.0).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_2).imag() == Approx(0.0).epsilon(1e-4));

    // Not sure what is wrong here. This doesn't work on appveyor, but works
    // on a 'real' Windows machine
    #ifndef _MSVC_LANG
    auto matrix_3 = lemon::Matrix3D(0, 1, 0, 0, 0, 1, 1, 0, 0);
    auto eigenvalues_3 = lemon::eigenvalues(matrix_3);
    CHECK(std::get<0>(eigenvalues_3).real() == Approx(-0.5).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_3).real() == Approx(-0.5).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_3).real() == Approx(1.0).epsilon(1e-4));
    CHECK(std::get<0>(eigenvalues_3).imag() == Approx(std::sqrt(3) / 2).epsilon(1e-4));
    CHECK(std::get<1>(eigenvalues_3).imag() == Approx(-std::sqrt(3) / 2).epsilon(1e-4));
    CHECK(std::get<2>(eigenvalues_3).imag() == Approx(0.0).epsilon(1e-4));
    #endif

    // Eigenvectors

    // Matrix 2 is actually the easiest case, so we start here

    auto eigenvectors_2 = lemon::eigenvectors(matrix_2, eigenvalues_2);
    CHECK(eigenvectors_2[0][0] == Approx(0.0 / std::sqrt(5.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[0][1] == Approx(2.0 / std::sqrt(5.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[0][2] == Approx(-1.0 / std::sqrt(5.0)).epsilon(1e-4));

    CHECK(eigenvectors_2[1][0] == Approx(-1.0 / std::sqrt(1.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[1][1] == Approx(0.0 / std::sqrt(1.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[1][2] == Approx(0.0 / std::sqrt(1.0)).epsilon(1e-4));

    CHECK(eigenvectors_2[2][0] == Approx(0.0 / std::sqrt(5.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[2][1] == Approx(1.0 / std::sqrt(5.0)).epsilon(1e-4));
    CHECK(eigenvectors_2[2][2] == Approx(2.0 / std::sqrt(5.0)).epsilon(1e-4));

    auto wikipedia_example = lemon::Matrix3D(3, 2, 6, 2, 2, 5, -2, -1, -4);
    auto eigenvalues_wkpd = lemon::eigenvalues(wikipedia_example);
    auto eigenvectors_wkpd = lemon::eigenvectors(wikipedia_example, eigenvalues_wkpd);

    CHECK(eigenvectors_wkpd[0][0] == Approx(-2.0 / std::sqrt(6.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[0][1] == Approx(-1.0 / std::sqrt(6.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[0][2] == Approx( 1.0 / std::sqrt(6.0)).epsilon(1e-4));

    CHECK(eigenvectors_wkpd[1][0] == Approx(-1.0 / std::sqrt(3.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[1][1] == Approx(-1.0 / std::sqrt(3.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[1][2] == Approx(1.0 / std::sqrt(3.0)).epsilon(1e-4));

    CHECK(eigenvectors_wkpd[2][0] == Approx(2.0 / std::sqrt(6.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[2][1] == Approx(1.0 / std::sqrt(6.0)).epsilon(1e-4));
    CHECK(eigenvectors_wkpd[2][2] == Approx(-1.0 / std::sqrt(6.0)).epsilon(1e-4));
}

TEST_CASE("Singular Value Decomposition") {

    auto matrix_2 = lemon::Matrix3D(2, 0, 0, 0, 3, 4, 0, 4, 9);
    auto svd = lemon::svd(matrix_2);

    CHECK(svd.U[0][0] == Approx(0.0).epsilon(1e-4));
    CHECK(-svd.U[0][1] == Approx(1.0).epsilon(1e-4));
    CHECK(svd.U[0][2] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.U[1][0] == Approx(0.447214).epsilon(1e-4));
    CHECK(svd.U[1][1] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.U[1][2] == Approx(0.894427).epsilon(1e-4));
    CHECK(svd.U[2][0] == Approx(0.894427).epsilon(1e-4));
    CHECK(svd.U[2][1] == Approx(0.0).epsilon(1e-4));
    CHECK(-svd.U[2][2] == Approx(0.447214).epsilon(1e-4));

    CHECK(svd.S[0][0] == Approx(11.0).epsilon(1e-4));
    CHECK(svd.S[0][1] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[0][2] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[1][0] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[1][1] == Approx(2.0).epsilon(1e-4));
    CHECK(svd.S[1][2] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[2][0] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[2][1] == Approx(0.0).epsilon(1e-4));
    CHECK(svd.S[2][2] == Approx(1.0).epsilon(1e-4));

    auto Vt = svd.V.transpose();
    CHECK(Vt[0][0] == Approx(0.0).epsilon(1e-4));
    CHECK(Vt[0][1] == Approx(0.447214).epsilon(1e-4));
    CHECK(Vt[0][2] == Approx(0.894427).epsilon(1e-4));
    CHECK(-Vt[1][0] == Approx(1.0).epsilon(1e-4));
    CHECK(Vt[1][1] == Approx(0.0).epsilon(1e-4));
    CHECK(Vt[1][2] == Approx(0.0).epsilon(1e-4));
    CHECK(Vt[2][0] == Approx(0.0).epsilon(1e-4));
    CHECK(Vt[2][1] == Approx(0.894427).epsilon(1e-4));
    CHECK(-Vt[2][2] == Approx(0.447214).epsilon(1e-4));

    auto m = svd.U * svd.S * svd.V.transpose();

    CHECK(m[0][0] == Approx(2.0).epsilon(1e-4));
    CHECK(m[0][1] == Approx(0.0).epsilon(1e-4));
    CHECK(m[0][2] == Approx(0.0).epsilon(1e-4));
    CHECK(m[1][0] == Approx(0.0).epsilon(1e-4));
    CHECK(m[1][1] == Approx(3.0).epsilon(1e-4));
    CHECK(m[1][2] == Approx(4.0).epsilon(1e-4));
    CHECK(m[2][0] == Approx(0.0).epsilon(1e-4));
    CHECK(m[2][1] == Approx(4.0).epsilon(1e-4));
    CHECK(m[2][2] == Approx(9.0).epsilon(1e-4));
}

TEST_CASE("Kabsch") {
    lemon::Coordinates in(100);
    lemon::Coordinates out(100);

    for (size_t row = 0; row < 100; ++row) {
        for (size_t col = 0; col < 3; ++col) {
            in[row][col] = std::log(2.0*static_cast<double>(row) + 10.0)/
                               std::sqrt(1.0*static_cast<double>(col) + 4.0) +
                               std::sqrt(static_cast<double>(col)*1.0)/
                               (static_cast<double>(row) + 1.0);
        }
    }

    auto s = lemon::Vector3D{ -5, 6, -27 };

    auto r = lemon::Matrix3D(
         -0.487179,   0.666667,   0.564103,
          0.871795,   0.333333,   0.358974,
          0.0512821,  0.666667,  -0.74359
    );

    for (size_t row = 0; row < 100; ++row) {
        out[row] = r*in[row];
        out[row] += s;
    }

    auto before_rmsd = lemon::rmsd(in, out);    

    auto affine = lemon::kabsch(in, out);

    CHECK(-affine.T[0] == Approx(5.0).epsilon(1e-4));
    CHECK(affine.T[1] == Approx(6.0).epsilon(1e-4));
    CHECK(-affine.T[2] == Approx(27.0).epsilon(1e-4));

    CHECK(-affine.R[0][0] == Approx(0.487179).epsilon(1e-4));
    CHECK(affine.R[0][1] == Approx(0.666667).epsilon(1e-4));
    CHECK(affine.R[0][2] == Approx(0.564103).epsilon(1e-4));
    CHECK(affine.R[1][0] == Approx(0.871795).epsilon(1e-4));
    CHECK(affine.R[1][1] == Approx(0.333333).epsilon(1e-4));
    CHECK(affine.R[1][2] == Approx(0.358974).epsilon(1e-4));
    CHECK(affine.R[2][0] == Approx(0.0512821).epsilon(1e-4));
    CHECK(affine.R[2][1] == Approx(0.666667).epsilon(1e-4));
    CHECK(-affine.R[2][2] == Approx(0.74359).epsilon(1e-4));

    lemon::align(out, affine);
    auto after_rmsd = lemon::rmsd(in, out);
    CHECK(after_rmsd < before_rmsd);
    CHECK(after_rmsd == Approx(0.0).epsilon(1e-6).epsilon(1e-4));
}

TEST_CASE("Error handling") {
    lemon::Coordinates a;
    CHECK(lemon::center(a) == lemon::Vector3D(0.0, 0.0, 0.0));

    lemon::Coordinates b;
    CHECK(lemon::rmsd(a, b) == 0);

    b.push_back({0.0, 0.0, 0.0});

    CHECK_THROWS_WITH(lemon::rmsd(a, b),
                      "lemon::rmsd(a,b): a.size() must equal b.size()");

    CHECK_THROWS_WITH(lemon::covariant(a, b),
                      "lemon::covariant(a,b): a.size() must equal b.size()");

    CHECK_THROWS_WITH(lemon::kabsch(a, b),
                      "lemon::kabsch(a,b): in.size() must equal out.size()");
}
