#ifndef LEMON_MATRIX_HPP
#define LEMON_MATRIX_HPP

#include <array>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles/types.hpp>
LEMON_EXTERNAL_FILE_POP

namespace lemon {

using Vector3D = chemfiles::Vector3D;
using Matrix3D = chemfiles::Matrix3D;
using Coordinates = std::vector<Vector3D>;

using chemfiles::cross;
using chemfiles::dot;

//! Error to issues with 
class matrix_error : public std::logic_error {
  public:
    explicit matrix_error(const std::string& what_arg)
        : std::logic_error(what_arg) {}
    explicit matrix_error(const char* what_arg)
        : std::logic_error(what_arg) {}
};

inline double trace(const Matrix3D& matrix) {
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}

template <typename Container = Coordinates>
inline Matrix3D covariant(const Container& a, const Container& b) {

    if (a.size() != b.size()) {
        throw matrix_error("lemon::covariant(a,b): a.size() must equal b.size()");
    }

    Matrix3D result = Matrix3D::zero();

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) { // transpose of b, not b
            for (size_t k = 0; k < a.size(); ++k) {    // We are stuck in 3 dimensions
                result[i][j] += a[k][i] * b[k][j]; // transpose of b, not b
            }
        }
    }

    return result;
}

template <typename Container = Coordinates>
inline double rmsd(const Container& a, const Container& b) {

    if (a.size() != b.size()) {
        throw matrix_error("lemon::rmsd(a,b): a.size() must equal b.size()");
    }

    if (a.size() == 0) {
        return 0.0;
    }

    double result = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        auto dist =  a[i] - b[i];
        result += dot(dist, dist);
    }
    return std::sqrt(result / static_cast<double>(a.size()));
}

template <typename Container = Coordinates>
inline Vector3D center(const Container& a) {
    Vector3D result = {0.0, 0.0, 0.0};

    if (a.size() == 0) {
        return result;
    }

    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i];
    }
    return result / static_cast<double>(a.size());
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
        // be greater than 1 as P is guaranteed to be either 0 or negative
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

    auto roots = cardano<Z>(a, b, c, d, eps);
    std::sort(roots.begin(), roots.end(), [](const Z& z1, const Z& z2){
        return std::abs(z1) < std::abs(z2);
    });

    return roots;
}

inline Vector3D cayley_hamilton(const Matrix3D& m, double e) {
    using std::fabs;
    auto col = (fabs(m[0][0]) > e || fabs(m[1][0]) > e || fabs(m[2][0]) > e)
            ? 0UL
            : 1UL;

    auto result = Vector3D{m[0][col], m[1][col], m[2][col]};
    return result / result.norm();
}

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

struct SingularValueDecomposition {
    Matrix3D U;
    Matrix3D S;
    Matrix3D V;
};

inline SingularValueDecomposition svd(const Matrix3D& m, double eps = 1e-10) {
    // clang-format off
    auto mt_m = m.transpose() * m;
    auto e_vals = eigenvalues(mt_m, eps);

    auto zero = std::complex<double>(0.0, 0.0);
    e_vals[0] = e_vals[0].real() > std::sqrt(eps) ? e_vals[0] : zero;
    e_vals[1] = e_vals[1].real() > std::sqrt(eps) ? e_vals[1] : zero;
    e_vals[2] = e_vals[2].real() > std::sqrt(eps) ? e_vals[2] : zero;

    auto rank = static_cast<int>(e_vals[0].real() > eps) +
                static_cast<int>(e_vals[1].real() > eps) +
                static_cast<int>(e_vals[2].real() > eps);

    auto e_vecs = eigenvectors(mt_m, e_vals, eps);

    // Create V with rows being the eigen vectors of m' * m
    // Note: the eigen vector function has already normalized these vectors
    auto Vt = Matrix3D(e_vecs[2][0], e_vecs[2][1], e_vecs[2][2],
                       e_vecs[1][0], e_vecs[1][1], e_vecs[1][2],
                       e_vecs[0][0], e_vecs[0][1], e_vecs[0][2]
    );
    auto V = Vt.transpose();

    Vector3D U_c0;
    Vector3D U_c1;
    Vector3D U_c2;

    if (rank == 2) {
        U_c0 = m * Vector3D{V[0][0], V[1][0], V[2][0]};
        U_c1 = m * Vector3D{V[0][1], V[1][1], V[2][1]};
        U_c2 = cross(U_c0, U_c1);
    } else {
        auto U = m * V;

        U_c0 = Vector3D{U[0][0], U[1][0], U[2][0]};
        U_c1 = Vector3D{U[0][1], U[1][1], U[2][1]};
        U_c2 = Vector3D{U[0][2], U[1][2], U[2][2]};
    }

    U_c0 = U_c0 / U_c0.norm();
    U_c1 = U_c1 / U_c1.norm();
    U_c2 = U_c2 / U_c2.norm();

    auto O = Matrix3D(U_c0[0], U_c1[0], U_c2[0],
                      U_c0[1], U_c1[1], U_c2[1],
                      U_c0[2], U_c1[2], U_c2[2]
    );

    auto S = Matrix3D(std::sqrt(e_vals[2].real()), 0.0, 0.0,
                      0.0, std::sqrt(e_vals[1].real()), 0.0,
                      0.0, 0.0, std::sqrt(e_vals[0].real())
    );

    return {O, S, V};
    // clang-format on
}

struct Affine {
    Vector3D T;
    Matrix3D R;
};

template <typename Container = Coordinates>
inline Affine kabsch(Container in, Container out, double eps = 1e-10) {

    if (in.size() != out.size()) {
        throw matrix_error("lemon::kabsch(a,b): in.size() must equal out.size()");
    }

    // Find the centroids then shift to the origin
    auto in_ctr = center(in);
    auto out_ctr = center(out);

    for (auto& i : in) {
        i -= in_ctr;
    }

    for (auto& i : out) {
        i -= out_ctr;
    }

    auto cov = covariant(out, in).transpose();
    auto SVD = svd(cov, eps);

    auto d = (SVD.V * SVD.U.transpose()).determinant();
    if (d > 0.0) {
        d = 1.0;
    } else {
        d = -1.0;
    }

    auto I = Matrix3D::unit();
    I[2][2] = d;
    auto R = SVD.V * I * SVD.U.transpose();

    return {std::move(out_ctr - R*in_ctr), std::move(R)};
}

template <typename Container = Coordinates>
inline void align(Container& in, const Affine& affine) {
    for (auto& c : in) {
        c -= affine.T;
        c = affine.R.transpose() * c;
    }
}

} // namespace lemon

#endif
