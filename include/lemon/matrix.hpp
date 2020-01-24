#ifndef LEMON_MATRIX_HPP
#define LEMON_MATRIX_HPP

#include <chemfiles/types.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

namespace lemon {

using Vector3D = chemfiles::Vector3D;
using Matrix3D = chemfiles::Matrix3D;

inline double trace(const Matrix3D& matrix) {
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}

inline Matrix3D covarience(const std::vector<Vector3D>& a,
                           const std::vector<Vector3D>& b) {

    Matrix3D result = Matrix3D::zero();

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) { // transpose of b, not b
            for (size_t k = 0; k < 3; ++k) {    // We are stuck in 3 dimensions
                result[i][j] += a[i][k] * b[j][k]; // transpose of b, not b
            }
        }
    }

    return result;
}

template <class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

template <typename Z = std::complex<double>>
inline std::array<Z, 3> cardano(double a, double b, double c, double d,
                                double eps = 1e-10) {
    
    auto a_3 = a * a * a;
    auto b_3 = b * b * b;

    // simplify into three new variables. see https://brilliant.org/wiki/cardano-method/
    auto Q = c / a - (b * b) / (a * a) / 3;

    auto R = (2 * b_3 / a_3 - 9 * b * c / (a * a)) / 27 + d / a;
    auto P = ((R * R) / 4) + ((Q * Q * Q) / 27);

    // all roots are the same and real
    if (std::fabs(Q) < eps && std::fabs(R) < eps && std::fabs(P) < eps) {
        auto root = Z(std::cbrt((d / a)), 0.0);
        return std::array<Z, 3>{{root, root, root}};
    }

    auto D = -(b / (3 * a));

    // The roots must be real as (R^2/4 - P) must be positive
    // However, this method may produce intermediate complex numbers, which
    // are avoided below via a cosine substitution
    if (P < eps) {
        auto i = std::sqrt(R * R / 4 - P);

        // For some reason, acos(1.0) and acos(-1.0) return nan instead of pi.
        // Otherwise, it is impossible for |R / (2 * sqrt(R^2 / 4 - P)) | to
        // be greater than 1 as P is gaurenteed to be either 0 or negative
        auto theta = std::acos(clamp(-R / (2 * i), -1.0 + eps, 1.0 - eps));

        auto m = std::cos(theta / 3);
        auto n = std::sqrt(3.0) * std::sin(theta / 3);

        auto j = std::cbrt(i);
        auto root_1 = Z(j * 2 * m + D, 0.0);
        auto root_2 = Z(-j * (m + n) + D, 0.0);
        auto root_3 = Z(-j * (m - n) + D, 0.0);

        return std::array<Z, 3>{{root_1, root_2, root_3}};
    }

    // P must be positive, so R^2 / 4 - P must be negative as P includes the
    // R^2 / 4 term. Therefore, we cannot use the cos(2*theta) trick as two of
    // the roots will be complex conjugates of each other. We calculate them
    // as shown on brilliant.org
    auto S = std::cbrt(-(R / 2) + std::sqrt(P));
    auto T = std::cbrt(-(R / 2) - std::sqrt(P));

    auto real_root = Z(D + (S + T), 0.0);

    auto complex_root = Z(D - (S + T) / 2, (S - T) * (std::sqrt(3.0) / 2));

    return std::array<Z, 3>{{real_root, complex_root, std::conj(complex_root)}};
}

template <typename Z = std::complex<double>>
inline std::array<Z, 3> eigenvalues(const Matrix3D& m, double eps = 1e-10) {
    // Setup the cubic polynomial
    auto a = 1.0;
    auto b = -trace(m);
    auto c = -0.5 * (trace(m * m) - trace(m) * trace(m));
    auto d = -m.determinant();

    return cardano<Z>(a, b, c, d, eps);
}

namespace {

Vector3D cayley_hamilton(const Matrix3D& m, double e) {
    using std::fabs;
    size_t col = (fabs(m[0][0]) > e || fabs(m[1][0]) > e || fabs(m[2][0]) > e)
            ? 0
            : 1;

    auto result = Vector3D{m[0][col], m[1][col], m[2][col]};
    return result / result.norm();
}

} // namespace

template <typename Z = std::complex<double>>
inline std::array<Vector3D, 3> eigenvectors(const Matrix3D& m,
                                            const std::array<Z, 3>& evals,
                                            double eps = 1e-10) {

    auto alpha1 = m - evals[0].real() * Matrix3D::unit();
    auto alpha2 = m - evals[1].real() * Matrix3D::unit();
    auto alpha3 = m - evals[2].real() * Matrix3D::unit();

    return std::array<Vector3D, 3>{cayley_hamilton(alpha2 * alpha3, eps),
                                   cayley_hamilton(alpha1 * alpha3, eps),
                                   cayley_hamilton(alpha1 * alpha2, eps)};
}

} // namespace lemon

#endif
